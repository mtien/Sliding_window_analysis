#!usr/bin/python

import sys, os, math

def reverseComplement(seq):
	complement=""	
	for i in seq:
		if( i =="A"):
			complement+="T"
		elif( i=="C"):
			complement+="G"
		elif(i=="G"):
			complement+="C"
		else:
			complement+="A"
	return complement[::-1]

def getGenomicSequence(loci_left, loci_right, loci_strand):
	file_genome= open(os.path.realpath("genome_files/CB15_SW/CP001340_2016_09_22.fna"))
	file_genome.readline()
	count=0
	seq=""
	
	length = loci_right - loci_left

	for line in file_genome:
		start= loci_left-count
		seq_line=line.strip()
		##find line where loci_left is located
		if( start <= len( seq_line ) and len(seq) <= 0 ):
			##start recording sequences till the end is reached			
			temp_seq=seq_line[start-1:]
			for i in range(len(temp_seq)):
				if( length-len(seq) > -1 ):
					seq+= temp_seq[i]

		elif( len(seq)>0 and length-len(seq) >-1 ):
			for i in range(len(seq_line)):
				if( length-len(seq) >-1 ):
					seq+= seq_line[i]	
		count+= len( seq_line )

	if(loci_strand=="+"):
		return seq
	else:
		return reverseComplement(seq)
		
def writeFASTA(headers, sequences, FASTA):
	f_out=open(FASTA, "w")
	for i in range(len(headers)):
		f_out.write(">"+ headers[i] + "\n")
		seq= sequences[i]
		temp_seq=""
		for base in seq:
			if(len(temp_seq) <70):
				temp_seq+=base
			else:
				f_out.write(temp_seq + "\n")
				temp_seq=base
		f_out.write(temp_seq+"\n")
	f_out.close()
