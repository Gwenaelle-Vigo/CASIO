#! /usr/bin/env python3

import sys

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


gff = open(sys.argv[1])
fasta_file = open(sys.argv[2])

### compute transcript size and store them in dictionary
transcriptSize={}

for name, seq in read_fasta(fasta_file):
	name=name.lstrip(">")
	transcriptSize["Parent="+name] = len(seq)


gene2Exclude = {}
gene2Keep = {}
# loop over the GFF from bedtools window
for li in gff:
	li = li.rstrip()
	# split the line according to the "tab"
	arrli = li.split("\t")
	gene1 = arrli[8].split(";")[0]
	gene2 = arrli[17].split(";")[0]
	if gene1 == gene2:
		continue
	## compare size
	sgene1 = transcriptSize[gene1]
	sgene2 = transcriptSize[gene2]
	if sgene1 > sgene2:
		## if longest transcript have already been exclude in previous comparison; then exclude also gene2 and pass
		if gene1 in gene2Exclude:
			gene2Exclude[gene2] = 1
			continue
		## if shorest gene is already un gene2Keep, then move it to gene2Exclude
		if gene2 in gene2Keep:
			del gene2Keep[gene2]
		gene2Keep[gene1] = 1
		gene2Exclude[gene2] = 1
	if sgene2 > sgene1:
		## if longest transcript have already been exclude in previous comparison; then exclude also gene2 and pass
		if gene2 in gene2Exclude:
			gene2Exclude[gene1] = 1
			continue
		## if shorest gene is already un gene2Keep, then move it to gene2Exclude
		if gene1 in gene2Keep:
			del gene2Keep[gene1]
		gene2Keep[gene2] = 1
		gene2Exclude[gene1] = 1
	## if gene has the same size
	if sgene1 == sgene2:
		if gene1 in gene2Exclude and gene2 in gene2Exclude:
			continue
		if gene1 not in gene2Keep and gene2 not in gene2Keep:
			gene2Keep[gene1] = 1
		## if both are in gene2Keep, one must excluded.
		if gene1 in gene2Keep and gene2 in gene2Keep:
			del gene2Keep[gene2]
			gene2Exclude[gene2] = 1


for gene in gene2Exclude:
	print(gene)

