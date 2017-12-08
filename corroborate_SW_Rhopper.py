#!usr/bin/python

import sys, os

rHopper=[]
file_tscript=open("result_files//NC_999999_transcripts_criterion.txt")
header=file_tscript.readline()
count= len(header.strip().split("\t"))
for line in file_tscript:
	info= line.strip().split("\t")
	if(len(info)>1):
		if( len(info) < count):
			info= [""]+info
		rHopper.append(info[6])
file_tscript.close()

count={key: 0 for key in rHopper}
file_in= open("data_files/POSvsNEG_p10_BMA1000_POS_windowInfo.txt")
header= file_in.readline()
##ID[0], BM[1], BMA[2], BMB[3], Fold[4], Log2[5], Pval[6], FDR[7], proximal_genes[8], left[9], right[10]
file_out= open("data_files/POSvsNEG_p10_BMA1000_POS_windowInfo_RhopperCongruent.txt", "w")
file_out.write(header)
file_out2= open("data_files/POSvsNEG_p10_BMA1000_POS_windowInfo_RhopperNotCongruent.txt", "w")
file_out2.write(header)
gene_not=["CCNA_R0069", "CCNA_R0066", "CCNA_R0065", "CCNA_R0087", "CCNA_R0084", "CCNA_R0083"]

for line in file_in:
	info=line.strip().split()
	if( info[8] != "" ):
		genes= info[8].split(",")[:-1]
		flag=False
		flag2= True
		for g in genes:
			if( g in rHopper ):
				flag=True
				count[g]+=1
			if( g in gene_not):
				flag2= False
		if(flag):
			file_out.write(line)
		elif(flag2):
			file_out2.write(line)

file_in.close()
file_out.close()
file_out2.close()

transcripts=[]
for r in rHopper:
	if(count[r] == 0):
		transcripts.append(r)

file_tscript=open("result_files/NC_999999_transcripts_criterion.txt")
header=file_tscript.readline()
file_out= open("data_files/NC_999999_transcripts_unique.txt", "w")
file_out.write(header)
for line in file_tscript:
	for t in transcripts:
		if(line.find(t) > 0):
			file_out.write(line)
file_tscript.close()
file_out.close()
