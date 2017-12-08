#!usr/bin/python

import sys, os
import FASTA_methods
	
if __name__ == "__main__":
	file_in=open("result_files/intarna_websrv_table_truncated.csv")
	pos_keys=range(1,113,1)
	gsrN_dic={key: 0 for key in pos_keys}
	gsrN="GCGAAGGCGGCGACGCCCCCCGACTTGAGACTGGCGAAAGCCAGTCCGCCGCCTTCGCACCCCTTTGTTGATGCGAACTCAGGCGGACCGGGTGACCGGACCGCCTTTCCTT"
	targets=[]
	sequences=[]
	
	targets_50=[]
	sequences_50=[]
	for line in file_in:
		info=line.split(";")
		
		gsrN_left= int(info[10])
		gsrN_right= int(info[11])
		
		count_range= range(gsrN_left, gsrN_right+1,1)
		for c in count_range:
			gsrN_dic[c]+=1
		
		third_line= info[-1].split("\\n")[3].split("_")
		seq=""
		for entry in third_line:
			if(entry!= ""):
				seq+=entry
		
		if(gsrN_right < 51):
			targets.append( info[0] )
			sequences.append(seq)
		else:
			targets_50.append( info[0] )
			sequences_50.append( seq)
	file_in.close()
	
	FASTA_methods.writeFASTA(targets, sequences, "result_files/IntaRNA_target_sequences_5prime.fasta")
	FASTA_methods.writeFASTA(targets_50, sequences_50, "result_files/IntaRNA_target_sequences_3prime.fasta")
	
	file_out=open("result_files/gsrN_binding_hotspots.txt", "w")
	file_out.write("position\tcount\tNuc\n")
	for p in pos_keys:
		file_out.write(str(p) +"\t"+ str(gsrN_dic[p])+"\t"+ gsrN[p-1] +"\n")
	file_out.close()
