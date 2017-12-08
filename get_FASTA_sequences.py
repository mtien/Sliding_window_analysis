#!usr/bin/python

import sys, os
import assign_window
import FASTA_methods

if __name__ == "__main__":
	dic= assign_window.make_gene_dic()
	
	file_in= open("result_files/Sliding_Window_analysis_results.txt")
	file_in.readline()
	
	##windows +-50
	window_dic={}
	
	gene_not=["CCNA_R0081"]
	for line in file_in:
		info=line.strip().split("\t")
		w_ID= info[0]
		left= int(info[1])
		right= int(info[2])
		strand= info[3]
		category= info[9]
		sense=info[10]
		gene_matches= info[8]
		min_left= int(info[11])
		max_right= int(info[12])
		g_strand=info[13]
		
		if( gene_matches not in gene_not):
			p_left= left-50
			p_right= right+50
			if(p_right-p_left > 2000):
				p_left= left+50
				p_right= right-50
			
			region= str(p_left) + "-" + str(p_right) + "(" + strand + "," + category+ "," + sense + ")"
			identity=gene_matches + ","+ str(left)
			if( sense == "S" ):
				window_dic[identity]= [p_left, p_right, strand, region]
			else:
				window_dic[identity+"A"]= [p_left, p_right, strand, region]
			
	file_in.close()
	
	##rHopper is all internal genes
	file_tscript=open("data_files/NC_999999_transcripts_unique.txt")
	header=file_tscript.readline()
	count= len(header.strip().split("\t"))
	for line in file_tscript:
		info= line.strip().split("\t")
		if(len(info)>1):
			if(len(info) < count):
				info= [""] + info
			
			syn=info[6]
			
			if( syn not in gene_not):
				left= int(dic[syn][0])
				right= int(dic[syn][1])
				strand= dic[syn][2]
			
				tracker=1
				while(right-left >2000):
					region= str(left) + "-" + str(left+1999) + "(" + strand + ")"
					window_dic[syn+","+str(tracker)]= [left, left+1999, strand, region]
					tracker+=1
					left+=1999
				region= str(left) + "-" + str(right) + "(" + strand + ")"
				window_dic[syn+","+str(tracker)]= [left, right, strand, region]
			
	file_tscript.close()
	
	windows= window_dic.keys()
	promo=""
	FASTA_headers=[]
	sequences=[]
	print len(windows)
	for w in windows:
		gene_info= window_dic[w]
		left=gene_info[0]
		right= gene_info[1]
		strand= gene_info[2]
		region= gene_info[3]
		promo=FASTA_methods.getGenomicSequence(left, right, strand)
		sequences.append(promo)
		
		header= w.split(",")[0] + "_" + str(left) + "_"+ str(right)
		FASTA_headers.append( header )
	FASTA_methods.writeFASTA(FASTA_headers, sequences, "result_files/IntaRNA_input.fasta")
	
