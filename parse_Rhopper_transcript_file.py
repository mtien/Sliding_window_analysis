#!usr/bin/python

import sys, os, math

def calculateLog2Fold(A, B):
	if(B >0):
		fold= A/B
	else:
		fold= A/1.0
	return math.log(fold,2)

if __name__ == "__main__":
	file_tscript=open("data_files/NC_999999_transcripts.txt")
	file_out= open("result_files/NC_999999_transcripts_criterion.txt", "w")
	header=file_tscript.readline()
	count= len(header.strip().split("\t"))
	header= header.strip()+ "\tLog2Fold\n"
	file_out.write(header)
	for line in file_tscript:
		info= line.strip().split("\t")
		if(len(info)>1):
			if( len(info) < count):
				info= [""]+info
			if( (float(info[27])<.05) and (float(info[17])>float(info[25])) and (float(info[17])>1000) ):
				log2Fold=calculateLog2Fold(float(info[17]), float(info[25]))
				if(  len(line.strip().split("\t")) < count):
					file_out.write("\t"+line.strip()+ "\t"+ str(log2Fold) +"\n")
				else:
					file_out.write(line.strip()+ "\t"+ str(log2Fold) +"\n")
	file_tscript.close()
	file_out.close()
