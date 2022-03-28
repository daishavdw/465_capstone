import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

infile = open(input_path)
outfile = open(output_path,'w')

for line in infile:
    line = line.strip()
    fields = line.split('\t')
    if fields[12] != 'NA':
        outfile.write(fields[0] + '\t' + fields[1] + '\t' +fields[2] + '\t' +  fields[3] + '\t' + fields[4] + '\t' + fields[7] +'\t' + fields[12] + '\n')
