"""

13-09-2023 KKS

This script takes input the metadata json file and creates a SAF file for featureCounts

"""

import json, sys

print(sys.argv)

with open(sys.argv[1],'r') as json_file:
	json_object=json.load(json_file)

print(json_object)


print("metadata.gene_name: %s" %(json_object['metadata.gene_name'][0]))
print("metadata.chromosome: %s" %(json_object['metadata.chromosome'][0]))

gene_name = json_object['metadata.gene_name'][0]

gene_chromosome = json_object['metadata.chromosome'][0]

gene_start = json_object['metadata.start'][0]

gene_end = json_object['metadata.end'][0]

gene_strand = ""
gene_anti_sense = ""
if json_object['metadata.strand'][0] > 0:
	gene_strand = "+"
	gene_anti_sense = "-"
else:
	gene_strand = "-"
	gene_anti_sense = "+"

print(f"gene_name: {gene_name} chromosome: {gene_chromosome} gene_start: {gene_start} gene_end: {gene_end} gene_strand: {gene_strand}")

# Create SAF file: the file name will be the gene name

filename=f"{gene_name}.saf"
print(filename)

SAF_file = open(filename, "w")

header = "GeneID\tChr\tStart\tEnd\tStrand\n"
SAF_file.write(header)

gene_feature = f"{gene_name}\t{gene_chromosome}\t{gene_start}\t{gene_end}\t{gene_strand}\n"
gene_feature_anti = f"{gene_name}anti\t{gene_chromosome}\t{gene_start}\t{gene_end}\t{gene_anti_sense}\n"
SAF_file.write(gene_feature)
SAF_file.write(gene_feature_anti)

SAF_file.close()