# 465_capstone

For filter_cadd.py, parse_sift.py, and pull_cadd_ranges.py: 
- Each of these scripts take in two arguements: input file and output file name

SIFT
- raw SIFT data can be found here: https://sift.bii.a-star.edu.sg/sift4g/public/Homo_sapiens/GRCh37.74/
- each of the files with a chromosome name + '.gz' were downloaded. These were merged using the following command: cat 1 2 3 .... > merged_sift
- This file was used as the input file for parse_sift.py which removed all of the rows that did not have a SIFT Score and removed all columns expect the chromosome, position, alternate allele, reference allele, gene name, and sift score.

CADD
- raw CADD data can be found here: https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz
- This data is the input for filter_cadd.py. This script returns a new file that removes rows with duplicate variant ID's but the same CADD PHRED score. It also removes all columns except chromosome, position, alternate allele, reference allele, gene name, and CADD PHRED score. 
- The output of filter_cadd.py is used as input for pull_cadd_ranges.py which returns the number of variants that fall within each range used in the histogram. This also returns the number of total variants in the file
- The number of unique genes present in the filtered file was determined by the following command:  awk '{print $6}' filtered_cadd | sort | uniq > unique_genes
