import sys

infilename = sys.argv[1]
infile = open(infilename, 'r')

gene_set = {}
output_lines = {}

header = infile.readline()
print(header.rstrip())

for line in infile:
	#print line
	cur_gene_id = line.split('\t')[0]
	#print line.split('\t')	
	cur_fpkm = float(line.split('\t')[9])
		
	if cur_gene_id not in gene_set:				
		gene_set[cur_gene_id] = cur_fpkm
		output_lines[cur_gene_id] = line
		
	else:
		existing_fpkm = gene_set[cur_gene_id]
		if cur_fpkm > existing_fpkm:		
			gene_set[cur_gene_id] = cur_fpkm
			output_lines[cur_gene_id] = line		

for k in sorted(output_lines):
	print(str(output_lines[k]).rstrip())
 		
infile.close()
