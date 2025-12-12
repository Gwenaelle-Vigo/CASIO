#! /usr/bin/env python3

import sys
import numpy as np


fp = open(sys.argv[1])
listCDS = open(sys.argv[2])

CDS2Ex = []
for line in listCDS:
	line = line.rstrip()
	CDS2Ex.append(line)



for line in fp:
	if line[0] == "#":
		continue
	line = line.rstrip()
	arrl = line.split("\t")
	geneName = arrl[8].split(";")[0]
	#geneName = geneName.lstrip("Parent=")
	#geneName = geneName.strip('"')
	#print(geneName)
	if geneName not in CDS2Ex:
		print(line)

