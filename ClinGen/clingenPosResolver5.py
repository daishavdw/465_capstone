import re #Used to search through data to pull required strings
import csv #Used to handle reading and writing of data
import requests, ssl #Used to peform REST API calls. The ssl package is required to make calls to servers over 'HTTPS'.
import time #Used to time requests to Clinvar server to prevent loss of data. Also used to calculate script runtimes.
import urllib.request, urllib.error #Used to download ClinGen database and inform of download failures. 
import os #Used for filepath checks and file deletion.
import sys #Used to handle commandline and end script at specific errors.
import math #Used for rounding in time calculation.

#Dictionary linking GRCh37 accession numbers to chromosome ids. 
GRCh37Chrom = {"NC_000001.10" : "1", "NC_000002.11" : "2", "NC_000003.11" : "3", 
                "NC_000004.11" : "4", "NC_000005.9" : "5", "NC_000006.11" : "6",
                "NC_000007.13" : "7", "NC_000008.10" : "8", "NC_000009.11" : "9",
                "NC_000010.10" : "10", "NC_000011.9" : "11", "NC_000012.11" : "12",
                "NC_000013.10" : "13", "NC_000014.8" : "14", "NC_000015.9" : "15",
                "NC_000016.9" : "16", "NC_000017.10" : "17", "NC_000018.9" : "18",
                "NC_000019.9" : "19", "NC_000020.10" : "20", "NC_000021.8" : "21",
                "NC_000022.10" : "22", "NC_000023.10" : "X", "NC_000024.9" : "Y",
                "NC_012920.1" : "MT"}


def findGRCh37Expression(expressions):
    key = "" #Refers to keys in the "GRCh37Chrom" dictionary created at beginning of file
    expression = "" #Refers to a single line string that contains info about the variation location.
    
    #Iterates through expressions list and GRCh37Chrom Dictionary and finds matches
    for e in expressions:
        for k in GRCh37Chrom:
            if k in e:
                key = k
                expression = e
    
    return([key, expression]) #Returns empty strings if there was not a GRCh37 Expression in the expressions cell of the row.

#Function for parsing VCF info from HGVS Expression
def parseFromExpression(key, expression):
    #Initialization of values for function use
    chrom = "null"
    pos = "null"
    ref = "null"
    alt = "null"
    parseMethod = "ERROR"

    #First, attempt to use Ensembl's Variant Recoder to discover the VCF info from the GRCh37 expression.

    baseURL = "http://grch37.rest.ensembl.org/variant_recoder/human/REPLACE/?vcf_string=1" #URL for the Ensembl GRCh37 REST server
    #Limit is 55,000 requests per hour so throttling shouldn't be an issue.
    #Server is generally slower than Clinvar and accepts batches of 200 in theory through a .POST call, but I haven't been able to get it to work.
    #.POST would be preferrable, but getting it to return VCF info through even the .GET function required some jury-rigging of the URL, going against API documentation.
    #Only to be used when ClinVar ID is unvailable or ClinVar does not have VCF Info.
    #Will fail to return VCF data if 'ref' or 'alt' alleles are larger than 5 kilobases.

    searchURL = baseURL.replace("REPLACE", expression) #Places the provided expression into correct position in URL

    attempts = 0
    maxAttempts = 5
    acquired = False
    while (acquired == False):
        try:
            variantRecoderRequest = requests.get(searchURL, timeout=5) #Performs REST request.
            acquired = True
        except requests.exceptions.Timeout:
            attempts+=1
            if (attempts == maxAttempts):
                print("CONNECTION FAILURE: Check status of http://grch37.rest.ensembl.org")
                sys.exit()

    recoded = re.findall(r'vcf_string:\s+.*-[ATGC]+', variantRecoderRequest.text) #Narrows down to VCF information
    if (len(recoded) == 1): #Checks to see if server returned valid VCF information
        vcf = re.findall(r'[\dYXMT]+-\d+-[ATGC]+-[ATGC]+', recoded[0])[0] #Narrows down to exact VCF info
        values = vcf.split('-') #Splits line containing VCF info into separate strings.
        values.append("RECODER")
        return(values)
    else: #If RECODER failed, will attempt to return flawed VCF info based on Expression alone. This can happen due to the exact positions of the variation being recorded as tenuous.
        #print(rows.index(r), key, expression)
        chrom = GRCh37Chrom[key] #Acquires chromosome number from GRCh37 Dictionary.

        if re.search(r'g\.\d+_', expression) != None:
            pos = "ISSUE" #Will be set as issue if position is more complex than 1 bp.
        elif re.search(r'\(\?', expression) != None:
            pos = "ISSUE" #Will be set as issue if expression indicates positions are unclear.
        else:
            pos = re.findall(r'g\.\d+', expression)[0].strip("g.") #Parses postion from expression

        genotype = expression[-3:] #Checks to see if mutation is a duplication or deletion
        if (genotype == "del") | (genotype == "dup"):
            alt = genotype
        elif genotype.find(">") != -1: #Checks to see if mutation is a substitution
            separate = genotype.split(">")
            ref = separate[0]
            alt = separate[1]
        else: #Raises issue if other mutation type.
            ref = "ISSUE"
            alt = "ISSUE"
    
        return([chrom, pos, ref, alt, "PARSED"])

def parseFromClinVar(idList):
    valuesDict = {}

    time.sleep(0.35) #Clinvar API restricts to 3 requests per second. Slows down script to prevent being locked out/loss of data.

    #These variables are required in this format for the 'requests' module to work
    fetchURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" #The base URL for performing an NCBI fetch.
    baseFetch = {"db" : "clinvar", "id" : "", "rettype" : "vcv", "is_variationid" : ""} #Dictionary containing modifiers for the fetch URL necessary for getting data from Clinvar.

    newFetch = baseFetch #Initializes a new fetch dictionary from the base.
    newFetch["id"] = idList #Modifys the new fetch dictionary to include the Clinvar ID from the file.

    attempts = 0
    maxAttempts = 5
    acquired = False
    while (acquired == False):
        try:
            response = requests.post(fetchURL, data = newFetch, timeout=10) #Performs the HTTP request.
            acquired = True
        except requests.exceptions.Timeout:
            attempts += 1
            if (attempts == maxAttempts):
                print("CONNECTION FAILURE: Check status of https://eutils.ncbi.nlm.nih.gov/")
                sys.exit()

    #print(response.text)
    #^Shows raw XML response from Clinvar

    #This block of code processes the raw response string from CLinVar so that it can be efficently parsed by regex
    rawResponseText = response.text
    singleLineResponseText = rawResponseText.replace("\n", " ").replace("\r", " ") #Removes all returns from response
    delimitedEntriesResponseText = singleLineResponseText.replace("</VariationArchive>", "\n\n") #Reintroduces returns only between separate entries

    entries = re.findall(r'VariationArchive.+', delimitedEntriesResponseText, re.U) #Separates entries into individual strings
    print("Entries found: " + str(len(entries)))

    for e in entries:
        
        #Narrows down string to position of ClinVar ID
        idPosition = re.findall(r'VariationArchive VariationID=\"\d+\"', e, re.U)
        if len(idPosition) == 1:
            id = re.findall(r'\d+', idPosition[0], re.U)[0] #Extracts exact ClinVar ID
        else:
            print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH") #Scream if ClinVar entry cannot be found by ID even if the id was searched for and returned.
        
        
        backReturn = e.replace("<", "\n<") #Separates the XML entries to be better parsed by regex
        #Narrows the response string down to the section containing information about the GRCh37 Assembly w/ VCF info
        clinvarGRCh37 = re.findall(r'<SequenceLocation Assembly=\"GRCh37\".+positionVCF[=\"\d \w+]+/>', backReturn, re.U)

        if (len(clinvarGRCh37) != 0): #Checks to see if Clinvar record does contain GRCh37 information

            chromEntry = re.findall(r'Chr=\"\w+\"', clinvarGRCh37[0], re.U)  #Narrows the string down to the section containing information about chromosome number.
            if (len(chromEntry) != 0): #Checks to see if Clinvar record does contain chromosome number information
                chrom = re.findall(r'\"[\dYyXxMmTt]+\"', chromEntry[0], re.U)[0].strip("\"") #Extracts exact chromosome number from string
                #This regular expression will pick up numbered chromosomes, Y chromosomes, X chromosomes, and mitochondrial chromosomes.

            posEntry = re.findall(r'positionVCF=\"\d+\"', clinvarGRCh37[0], re.U) #Narrows the string down to the section containing information about position on chromosome.
            #Currently will exclude mutations longer than 1 bp
            if (len(posEntry) != 0): #Checks to see if Clinvar record does contain information about position on chromosome.
                pos = re.findall(r'\d+', posEntry[0], re.U)[0] #Extracts exact position from string

            refEntry = re.findall(r'referenceAlleleVCF=\"[ATCG]+\"', clinvarGRCh37[0], re.U) #Narrows the string down to the section containing reference allele information.
            #Regex needs work to expand mutation diversity.
            if (len(refEntry) != 0): #Checks to see if Clinvar record does contain reference allele information.
                ref = re.findall(r'\"[ATCG]+\"', refEntry[0], re.U)[0].strip("\"") #Extracts exact reference allele from string.

            altEntry = re.findall(r'alternateAlleleVCF=\"[ATCG]+\"', clinvarGRCh37[0], re.U) #Narrows the string down to the section containing alternate allele information.
            if (len(altEntry) != 0): #Checks to see if Clinvar record does contain alternate allele information.
                    alt = re.findall(r'\"[ATCG]+\"', altEntry[0], re.U)[0].strip("\"") #Extracts exact alternate allele from string.
            
            valuesDict[id] = [chrom, pos, ref, alt, "CLINVAR"]         
        
        #If ClinVar does not provide the VCF information for the variant, send back with values indicating parsing required
        else:
            #Default to null values
            valuesDict[id] = ["null", "null", "null", "null", "PARSING REQUIRED"]

    #Searches through provided ID list to make sure every ID returned an entry. If not, flag for parsing.
    #This appears to happen rarely were a provided Clinvar ID returns no data from the API whatsoever.
    for x in idList:
        if x not in valuesDict:
            valuesDict[x] = ["null", "null", "null", "null", "PARSING REQUIRED"]


    return(valuesDict)


#------------BEGIN CMD LINE ARGUMENT HANDLING------------

if len(sys.argv) > 2:
    filePath = sys.argv[2] #Extracts new file path from second command line argument.
else:
    filePath = "clinGenPathogenicity.csv" #Sets default file path if none is specified in args

newDownload = -1 #Initialization. Controls whether new data is downloaded. 1 = Yes, 0 = No.
if len(sys.argv) == 1: #Checks to see if someone has input the cmd line argument for whether new data should be downloaded.
    while newDownload == -1: #Asks question until valid input is given.
        answer = input("Would you like to download fresh data? (Y/N)\n")

        #Parses terminal answer and sets 'newDownload' value accordingly
        if (re.search(r'[Yy]', answer) != None):
            newDownload = 1
        elif (re.search(r'[Nn]', answer) != None):
            newDownload = 0
else:
    newDownload = int(sys.argv[1]) #Gets 'newDownload' value from cmd line argument if given.


#------------END CMD LINE ARGUMENT HANDLING------------


#------------BEGIN FILE DOWNLOAD HANDLING-------------
if ((newDownload == 1) | (os.path.exists(filePath) == False)): #Determines if new download was requested or if file doesn't currently exists under expected name.
    if os.path.exists(filePath): #Determines if old version exists in directory
        try:
            urllib.request.urlopen("https://erepo.clinicalgenome.org/evrepo/api/classifications/all") #Attempts to connect to URL where file is known to be.
        except urllib.error.HTTPError: # Raises error if it cannot connect to URL. File has likely moved to somewhere else on their website.
            print("ERROR: Cannot find file at specified URL. Check ClinGen website for potentially updated URL. Old data will not be deleted.")
            sys.exit()
        os.remove(filePath) #Removes old file to prevent duplicate data. Should not occur if new file could not be found at specified URL.
    elif newDownload == 0: #Informs user that new data will have to be downloaded if they did not request it.
        print("ERROR: No existing data found. Downloading fresh data.")

    try:
        downloadStartTime = time.time()
        print("Download Start: " + time.ctime(downloadStartTime))
        
        print("Downloading...")
        
        urllib.request.urlretrieve("https://erepo.clinicalgenome.org/evrepo/api/classifications/all", filePath) #Initiates download.
        #Server is slow and will take a few minutes. New file is created instantly and data is saved to it as it is downloaded.
        
        downloadEndTime = time.time()
        downloadDuration = downloadEndTime - downloadStartTime
        downDurationMinutes = math.floor(downloadDuration / 60)
        downDurationSeconds = round(downloadDuration % 60)
        print("Download Complete: " + time.ctime(downloadEndTime) + " - " + str(downDurationMinutes) + ":" + str(downDurationSeconds) + " elapsed.")
    except urllib.error.HTTPError: #Raises error if it cannot connect to URL. File has likely moved to somewhere else on their website.
        print("ERROR: Cannot find file at specified URL. Check ClinGen website for potentially updated URL.")
        sys.exit()


#------------END FILE DOWNLOAD HANDLING-------------

#------------BEGIN POSITION RESOLUTION-------------
resolutionStartTime = time.time()
print("Position Resolution Begin: " + time.ctime(resolutionStartTime))
rows = [] #Initialization of list to hold file rows.

#Opens stored file.
with open(filePath, encoding='utf-8') as clinGen: #"encoding='utf-8'" flag must be present. For whatever reason, file cannot be read without it.
    pathoReader = csv.reader(clinGen, delimiter = "\t")

    #Goes through each row in file, parsing and adding it to the list of rows.
    for row in pathoReader:
        rows.append(row)

#Inserts new columns into the header row. Inserts occur in reverse order, each at position 1, to account for the shifting of columns to the right after insert.
rows[0].insert(1, "parseMethod")
rows[0].insert(1, "alt")
rows[0].insert(1, "ref")
rows[0].insert(1, "pos")
rows[0].insert(1, "chrom")

existingRows = []

#Opens exisiting file
with open("clinGenPathogenicityFixedCoordinates.csv", 'r', newline="\n", encoding='utf-8') as existingFile:
    fileReader = csv.reader(existingFile)

    
    for row in fileReader:
        existingRows.append(row)

    existingFile.close

allClinVarIds = [] #Will contain a list of all ClinVar IDs in file.

if (len(existingRows) == len(rows)):
    print("No new rows found. Resolution complete.")
else:
    print(str(len(rows) - len(existingRows)) + " new rows found.")
    for r in rows[1:]:
        #Initializes new values. Values for specific rows that could not be found can be identified by this initial value.
        chrom = "null"
        pos = "null"
        ref = "null"
        alt = "null"
        parseMethod = "ERROR"
    
        #The ClinGen file will include multiple expressions of different formats for each listed variation. Often including GRCh37 and GRCh38. These are in the same cell and thus the GRCh37 expression must be searched for.
        # Splits the string in the cell containing expressions into multiple strings, each containing a single expression.
        expressions = r[3].split(", ")
        # Finds the GRCh37 Expression, if it exists. If not, returns empty strings.
        results = findGRCh37Expression(expressions)
        key = results[0]
        expression = results[1]
    
        clinvarID = r[1]  # Assigns ClinVar ID from proper cell.
        variationID = r[0]
        expressionEntry = r[3]

        addNewData = True


        if (rows.index(r) < len(existingRows)):
            if (existingRows[rows.index(r)][0] == variationID) & (existingRows[rows.index(r)][6] == clinvarID) & (existingRows[rows.index(r)][8] == expressionEntry):
                r.insert(1, existingRows[rows.index(r)][5])
                r.insert(1, existingRows[rows.index(r)][4])
                r.insert(1, existingRows[rows.index(r)][3])
                r.insert(1, existingRows[rows.index(r)][2])
                r.insert(1, existingRows[rows.index(r)][1])
                addNewData = False
        
        if (addNewData == True):
            #First attempt to obtain data using ClinVar API
            if clinvarID != '-':  # Makes sure that row has a valid CLinvar ID
    
                allClinVarIds.append(clinvarID)  # Adds Clinvar ID to master list.
    
            #If no ClinVar ID is given, attempt will be made to parse from expression, given that a valid expression exists.
            # If the GRCh37 expression was found, extract information from it.
            elif (key != "") & (expression != ""):
                #Get values using expression parser function
                values = parseFromExpression(key, expression)
                chrom = values[0]
                pos = values[1]
                ref = values[2]
                alt = values[3]
                parseMethod = values[4]
    
                #Inserts new values into the row. Inserts occur in reverse order, each at position 1, to account for the shifting of columns to the right after insert.
                r.insert(1, parseMethod)
                r.insert(1, alt)
                r.insert(1, ref)
                r.insert(1, pos)
                r.insert(1, chrom)
    
                print(rows.index(r), chrom, pos, ref, alt, parseMethod)
            else:
                #Default to null values
                #Inserts new values into the row. Inserts occur in reverse order, each at position 1, to account for the shifting of columns to the right after insert.
                r.insert(1, parseMethod)
                r.insert(1, alt)
                r.insert(1, ref)
                r.insert(1, pos)
                r.insert(1, chrom)
                print(rows.index(r), chrom, pos, ref, alt, "ERROR")
    
    remainingIDs = len(allClinVarIds)
    #print(remainingIDs)
    clinIDList = []  # List to be used to stage correct amount of clinvar ids for an API call
    # {ClinVar ID : [chrom, pos, ref, alt, parseMethod]} Will be used to store combined results of Clinvar API calls.
    allValuesDict = {}
    for i in allClinVarIds:  # Iterates through all Clinvar IDs
        clinIDList.append(i)  # Adds ID to staging list
    
        #ClinVar's efetch function takes a maximum of 200 IDs
        if (len(clinIDList) == 199):  # Makes sure staging list does not get larger than 200
            # Performs API call, returns segment of 'allValuesDict'
            valuesDict = parseFromClinVar(clinIDList)
    
            #Appends Clinvar API results to combinded dictionary
            for d in valuesDict:
                allValuesDict[d] = valuesDict[d]
    
            clinIDList.clear()  # Clears the staging list to start again
            remainingIDs = remainingIDs - 199
            print("IDs remaining: " + str(remainingIDs))
    
    #Perform final ClinVar API call using remaining IDs
    valuesDict = parseFromClinVar(clinIDList)
    
    #Adds remaining results to combined dictionary
    for d in valuesDict:
        allValuesDict[d] = valuesDict[d]
    
    clinIDList.clear()
    remainingIDs = 0
    print("CLINVAR API CALLS DONE")
    nullValues = ["DictMatchError", "DictMatchError",
                "DictMatchError", "DictMatchError", "DictMatchError"]
    
    #Merges data from ClinVar API dictionary with row data.
    for r in rows[1:]:  # Interates through rows.
        # Makes sure to use rows that have not yet had VCF info added (Rows containing Clinvar IDs).
        if len(r) < len(rows[0]):
            clinvarID = r[1]  # Assigns ClinVar ID from proper cell.
            values = nullValues
            #Iterates through combined dictionary and matches API data with row data.
            for id in allValuesDict:
                if id == clinvarID:
                    values = allValuesDict[id]
                    break
    
            # Determines if API call failed to get the data for this specific entry.
            if values[4] == "PARSING REQUIRED":
                #Extracts GRCh37 expression for parsing.
                expressions = r[3].split(", ")
                results = findGRCh37Expression(expressions)
                key = results[0]
                expression = results[1]
    
                if (key == "") | (expression == ""):
                    values[0] = "null"
                    values[1] = "null"
                    values[2] = "null"
                    values[3] = "null"
                    values[4] = "ERROR"
                else:
                    # Performs parse using Ensembl Variant Recoder API
                    values = parseFromExpression(key, expression)
            elif values != nullValues:
                # Pulls VCF info from combined dictionary if ClinVAr API call was successful.
                values = allValuesDict[id]
    
            #Assigns values from results.
            chrom = values[0]
            pos = values[1]
            ref = values[2]
            alt = values[3]
            parseMethod = values[4]
    
            #Inserts new values into the row. Inserts occur in reverse order, each at position 1, to account for the shifting of columns to the right after insert.
            r.insert(1, parseMethod)
            r.insert(1, alt)
            r.insert(1, ref)
            r.insert(1, pos)
            r.insert(1, chrom)
    
            print(rows.index(r), chrom, pos, ref, alt, parseMethod)
    
    manualRows = []
    # "encoding='utf-8'" flag must be present. For whatever reason, file cannot be read without it.
    with open("clinGenPathogenicity_ManualFix.csv", encoding='utf-8') as manual:
        fixer = csv.reader(manual)
    
        for f in fixer:
            manualRows.append(f)
    
    
    for r in rows[1:]:
        if ((r[5] == "PARSED") | (r[5] == "ERROR")):
            for m in manualRows[1:]:
                if (m[1] == r[6]) & (m[2] == r[0]) & (m[3] == r[8]):
                    r[1] = m[4]
                    r[2] = m[5]
                    r[3] = m[6].upper()
                    r[4] = m[7].upper()
                    r[5] = "MANUAL"
   

    #Opens file under same name and writes all lines (old+new)
    with open("clinGenPathogenicityFixedCoordinates.csv", 'w', newline="\n", encoding='utf-8') as newFile:
        fileWriter = csv.writer(newFile)

        for r in rows:
            fileWriter.writerow(r)

        newFile.close
    
    manual.close


#Calculates total computation time.
resolutionEndTime = time.time()
resolutionDuration = resolutionEndTime - resolutionStartTime
resDurationMinutes = math.floor(resolutionDuration / 60)
resDurationSeconds = round(resolutionDuration % 60)
print("Position Resolution Complete: " + time.ctime(resolutionEndTime) + " - " + str(resDurationMinutes) + ":" + str(resDurationSeconds) + " elapsed.")
clinGen.close
