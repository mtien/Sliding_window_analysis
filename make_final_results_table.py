#!usr/bin/python

import sys, os
import assign_window

if __name__ == "__main__":
	counter=0
	dic= assign_window.make_gene_dic()
	## syn: [gene, region, product]
	file_in= open("result_files/NC_999999_transcripts_criterion.txt")
	header=file_in.readline()
	count= len(header.strip().split("\t"))
	
	out_dic={}
	for line in file_in:
		info= line.strip().split("\t")
		if(len(info)>1):
			if( len(info) < count):
				info= [""]+info
			syn= info[6]
			info2= dic[syn]
			region= info2[0] + "-" + info2[1] + " (" + info2[2] + ")"
			##syn: [gene_name, enrichment_value, identification_method, region(s), product]
			out_dic[syn]= [info2[4], info[-1][:4], "Rockhopper", region, info2[5]]
			counter+=1
	file_in.close()
	
	file_in= open("result_files/Sliding_Window_analysis_results.txt")
	file_in.readline()
	for line in file_in:
		info=line.strip().split()
		left= info[1]
		right=info[2]
		strand=info[3]
		catag= info[9]
		sense= info[10]
		region= left+ "-" + right + " (" + strand + ", " + catag+ ", " + sense + ")"
		enrich= info[6][:4]
		if( info[8].find("CCNA") > -1):
			genes= info[8].split(",")
			if(len(genes) > 1):
				genes= genes[:-1]
			if(len(genes)==1):
				g=genes[0]
				if(g in out_dic.keys()):
					out_dic[g][1]+=", " + enrich
					out_dic[g][2]+=", Sliding window"
					out_dic[g][3]+=", " + region
				else:
					name= dic[g][4]
					prod= dic[g][5]
					out_dic[g]= [name, enrich, "Sliding window", region, prod]
				counter+=1
			else:
				g= genes[0]+ ","+ genes[1]
				name= dic[genes[0]][4] + ","+ dic[genes[1]][4]
				prod= dic[genes[0]][5] + ", "+ dic[genes[1]][5]
				out_dic[g]= [name, enrich, "Sliding window", region, prod]
				counter+=1
	file_in.close()
	
	print counter
	##file_out= open("result_files/Final_result_table.txt", "w")
	##file_out.write("Gene Locus ID\tGene Name\tEnrichment Value (Log2Fold)\t Identification Method\tRegions(s)\tDescription\n")
	check_genes= out_dic.keys()
	for c in check_genes:
		line= c
		for entry in out_dic[c]:
			line+="\t" + entry
		line+="\n"
		##file_out.write(line)
	##file_out.close()
