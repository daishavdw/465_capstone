import sys

infile_path = sys.argv[1]
outfile_path = sys.argv[2]

infile = open(infile_path)
outfile = open(outfile_path,'w')


less_than_10 = 0
less_than_20 = 0
less_than_30 = 0
less_than_40 = 0
less_than_50 = 0
less_than_60 = 0
less_than_70 = 0
less_than_80 = 0
less_than_90 = 0
less_than_100 = 0
count = 0
first = True
for line in infile:
    if first:
        first = False
    else:
        count += 1
        line = line.strip()
        fields = line.split('\t')
        if float(fields[6]) < 10:
            less_than_10 += 1
        elif float(fields[6]) < 20:
            less_than_20 += 1
        elif float(fields[6]) < 30:
            less_than_30 += 1
        elif float(fields[6]) < 40:
            less_than_40 += 1
        elif float(fields[6]) < 50:
            less_than_50 += 1
        elif float(fields[6]) < 60:
            less_than_60 += 1 
        elif float(fields[6]) < 70:
            less_than_70 += 1
        elif float(fields[6]) < 80:
            less_than_80 += 1
        elif float(fields[6]) < 90:
            less_than_90 += 1
        elif float(fields[6]) < 100:
            less_than_100 += 1  
    
outfile.write("less than 10: " + str(less_than_10) + '\n')
outfile.write("less than 20: " + str(less_than_20) + '\n')
outfile.write("less than 30: " + str(less_than_30) + '\n')
outfile.write("less than 40: " + str(less_than_40) + '\n')
outfile.write("less than 50: " + str(less_than_50) + '\n')
outfile.write("less than 60: " + str(less_than_60) + '\n')
outfile.write("less than 70: " + str(less_than_70) + '\n')
outfile.write("less than 80: " + str(less_than_80) + '\n')
outfile.write("less than 90: " + str(less_than_90) + '\n')
outfile.write("less than 100: " + str(less_than_100) + '\n')
outfile.write("count: " + str(count) + '\n')
