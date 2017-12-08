#!usr/bin/python

import sys, os
import numpy

file_in= open("data_files/RPKM_compiled_slidingWindow.txt")
header=file_in.readline()
file_out= open("data_files/RPKM_compiled_slidingWindow_removeHighVariantWindows.txt", "w")
file_out.write(header)
for line in file_in:
	info=line.strip().split("\t")
	gene_ID= info[0]
	POS_RNA=[]
	POS_2=[]
	for i in range(6,9):
		POS_RNA.append(int(info[i]))
	for i in range(7,9):
		POS_2.append(int(info[i]))
	check=numpy.mean(POS_2)
	if( check <1.0):
		check = 1.0
	if( numpy.mean(POS_RNA)/check < 1.5):
		file_out.write( line)
	
file_in.close()
file_out.close()
