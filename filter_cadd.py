
import  gzip
from io import TextIOWrapper
from gzip import GzipFile
import sys


infile_path = sys.argv[1]
outfile_path = sys.argv[2]

infile = open(infile_path, 'rb')
gzipped = GzipFile(None, 'rb', fileobj=infile)
data = TextIOWrapper(gzipped)

outfile = open(outfile_path,'w')

old_variant_id = ""
outfile.write("variant_id\tchrom\tpos\tref\talt\tgene_name\tPHRED\n")
for line in data:
	if line[0] != "#":
		line = line.strip()
		fields = line.split('\t')
		variant_id = fields[0] + ":" + fields[1] + ":" + fields[2] + ":" + fields[3]
		if variant_id != old_variant_id:
			old_variant_id = variant_id
			outfile.write(variant_id + '\t' + fields[0] + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[20] +  '\t' + fields[115])
			outfile.write('\n')	

