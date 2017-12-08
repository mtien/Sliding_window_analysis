#!usr/bin/python

import sys, os, math

def match_window(g_dic, window_info):
	left=int(window_info[0])
	right=int(window_info[1])
	genes= g_dic.keys()
	proximity=[]
	for g in genes:
		left_g = int(g_dic[g][0])
		right_g = int(g_dic[g][1])
		if( (left_g-200 < left) and (right_g+200 > right) ):
			proximity.append(g)
	return proximity

def find_window(window):
	file_window_name="genome_files/CB15_SW/CP001340_2016_09_22_SW_25.ptt"
	if(window=="CCNA_R0081"):
		file_window_name="genome_files/CB15_SW/CP001340_2016_09_22_SW_25.rnt"
	file_window=open(file_window_name)
	file_window.readline()
	file_window.readline()
	file_window.readline()
	window_info=["0","0","0"]
	for line in file_window:
		info= line.strip().split()
		synonym=info[5]
		if(window ==synonym ):
			info= line.strip().split()
			location=info[0]
			left=location.split("..")[0]
			right=location.split("..")[1]
			strand=info[1]
			length=info[2]
			PID=info[3]
			gene=info[4]
			Code=info[6]
			COG=info[7]
			Product=info[8]
			window_info=[left,right,synonym]
			break
	file_window.close()
	return window_info

def make_gene_dic(PPT= "genome_files/CP001340_genome_2016_09_22/CP001340_2016_09_22.ptt", RNT= "genome_files/CP001340_genome_2016_09_22/CP001340_2016_09_22.rnt"):
	
	gene_dic={}
	
	file_ppt= open(os.path.realpath(PPT))
	file_ppt.readline()
	file_ppt.readline()
	file_ppt.readline()
	for line in file_ppt:
		info= line.strip().split("\t")
		location=info[0]
		left=location.split("..")[0]
		right=location.split("..")[1]
		strand=info[1]
		length=info[2]
		PID=info[3]
		gene=info[4]
		synonym=info[5]
		Code=info[6]
		COG=info[7]
		product=info[8]
		gene_dic[synonym]=[left,right,strand,length,gene,product]
	file_ppt.close()
	
	file_rnt= open(os.path.realpath(RNT))
	file_rnt.readline()
	file_rnt.readline()
	file_rnt.readline()
	for line in file_rnt:
		info= line.strip().split("\t")
		location=info[0]
		left=location.split("..")[0]
		right=location.split("..")[1]
		strand=info[1]
		length=info[2]
		PID=info[3]
		gene=info[4]
		synonym=info[5]
		Code=info[6]
		COG=info[7]
		product=info[8]
		gene_dic[synonym]=[left,right,strand,length,gene,product]
	file_rnt.close()
	
	return gene_dic

def determine_position(strand, g_left, g_right, w_left, w_right):
	if(strand=="+"):
		##upstream
		if(w_right < g_left ):
			return "U", abs(g_left-w_right)
		##downstream
		elif( w_left > g_right):
			return "D", abs(w_left-g_right)
		##internal
		else:
			return "I", w_left-g_left
	else:
		##upstream
		if(g_right < w_left ):
			return "U", abs(w_left-g_right)
		##downstream
		elif( g_left > w_right):
			return "D", abs(g_left-w_right)
		##internal
		else:
			return "I", g_right-w_right
		
def combineWindows(geneID_dictionary):
	sorted_windows= sorted(geneID_dictionary)
	new_dic= {}
	
	working_window=sorted_windows[0]
	working_left= geneID_dictionary[working_window][0]
	working_right= geneID_dictionary[working_window][1]
	
	previous_window=sorted_windows[0]
	for s in range(1,len(sorted_windows)):
		window=sorted_windows[s]
		if( window-previous_window == 1):
			working_right= geneID_dictionary[window][1]
		else:
			window_ID= "W_" + str(working_window)
			new_dic[window_ID]=[working_left,working_right]
			working_window=window
			working_left= geneID_dictionary[window][0]
			working_right= geneID_dictionary[window][1]
		previous_window= window
	window_ID= "W_" + str(working_window)
	new_dic[window_ID]=[working_left,working_right]
	
	return new_dic

def resolveMultipleInternalGenes(position_dictionary):
	pos_keys= position_dictionary.keys()
	pos=0
	neg=0
	left_most_gene=pos_keys[0]
	right_most_gene=pos_keys[0]
	min_left=position_dictionary[left_most_gene][3]
	max_right=position_dictionary[left_most_gene][4]
	for p in pos_keys:
		entries= position_dictionary[p]
		###position, magnitude, strand, left, right
		p_left= entries[3]
		p_right=entries[4]
		p_strand= entries[2]
		if(p_strand=="+"):
			pos+=1
		else:
			neg+=1
		if(p_left < min_left):
			min_left=p_left
			left_most_gene=p
		if(p_right > max_right):
			max_right=p_right
			right_most_gene=p
	if(pos>neg):
		return left_most_gene, "+"
	else:
		return right_most_gene, "-"

def match_window_get_strand(g_dic, left, right):
	genes= g_dic.keys()
	proximity=[]
	for g in genes:
		left_g = int(g_dic[g][0])
		right_g = int(g_dic[g][1])
		if( right-left < right_g-left_g ):
			##window is smaller than gene
			if( (left_g-200 < left) and (right_g+200 > right) and (g not in proximity)):
				proximity.append(g)
		else:
			##window is bigger than gene
			##check internal
			if( left< left_g and right_g < right and (g not in proximity)):
				proximity.append(g)
			elif( (left_g-200 < left) and (right_g+200 > right) and (g not in proximity) ):
				proximity.append(g)
			
	
	if(len(proximity) > 1):
		##several proximity genes are detected, make a position dictionary
		pos_dic={}
		for p in proximity:
			##left, right, strand, length, gene, product
			entries= g_dic[p]
			g_strand= entries[2]
			left_g= int(entries[0])
			right_g= int(entries[1])
			position, magnitude= determine_position(g_strand, left_g, right_g, left, right)
			pos_dic[p]=[position, magnitude, g_strand, left_g, right_g]
		
		gene_internal=""
		internal_count= 0
		pd= pos_dic.keys()
		for p in pd:
			if(pos_dic[p][0] == "I"):
				gene_internal= p
				internal_count+=1
		
		if(internal_count == 1):
			return gene_internal, pos_dic[gene_internal][2]
		elif(internal_count < 1):
			min_mag=250
			for p in pd:
				mag= pos_dic[p][1]
				if(mag < min_mag):
					min_mag = mag
					gene_internal=p
			return gene_internal, pos_dic[gene_internal][2]
		else:
			## TOO MANY INTERNAL GENES!!!
			return resolveMultipleInternalGenes(pos_dic)
			
	elif( len(proximity) == 1 ):
		gene_ID= proximity[0]
		##left, right, strand, length, gene, product
		entries= g_dic[gene_ID]
		g_strand= entries[2]
		left_g= int(entries[0])
		right_g= int(entries[1])
		position, magnitude= determine_position(g_strand, left_g, right_g, left, right)
		return gene_ID, g_strand
	else:
		## NO GENES in proximity
		return "NA", "+"
	
def combineWindowsDic(geneID_dictionary):
	sorted_windows= sorted(geneID_dictionary)
	new_dic= {}
	
	working_window=sorted_windows[0]
	working_left= geneID_dictionary[working_window][0]
	working_right= geneID_dictionary[working_window][1]
	working_meanA= geneID_dictionary[working_window][2]
	working_meanB= geneID_dictionary[working_window][3]
	
	previous_window=sorted_windows[0]
	window_log=""
	window_count=1
	for s in range(1,len(sorted_windows)):
		window=sorted_windows[s]
		if( window-previous_window == 1):
			working_right= geneID_dictionary[window][1]
			working_meanA+= geneID_dictionary[window][2]
			working_meanB+= geneID_dictionary[window][3]
			window_log+="W_" + str(window) +","
			window_count+=1
		else:
			window_ID= "W_" + str(working_window)
			new_dic[window_ID]=[working_left,working_right,working_meanA, working_meanB, window_log, window_count]
			working_window=window
			working_left= geneID_dictionary[window][0]
			working_right= geneID_dictionary[window][1]
			working_meanA= geneID_dictionary[window][2]
			working_meanB= geneID_dictionary[window][3]
			window_log=""
			window_count=1
		previous_window= window
	window_ID= "W_" + str(working_window)
	new_dic[window_ID]=[working_left,working_right, working_meanA, working_meanB,window_log,window_count]
	
	return new_dic
	
def calculateLog2Fold(A, B):
	fold=1.0
	if(B >0):
		fold= float(A)/float(B)
	else:
		fold= float(A)/1.0
	return math.log(fold,2)

def getCombinedWindowInfo(geneID_dic):
	file_ptt= open(os.path.realpath("genome_files/CB15_SW/CP001340_2016_09_22_SW_25.ptt") )
	species = file_ptt.readline()
	ptt_Count= file_ptt.readline()
	header = file_ptt.readline()
	geneIDs= geneID_dic.keys()
	
	window_Dic= {}
	for line in file_ptt:
		for gene in geneIDs:
			info=line.strip().split("\t")
			check=info[5]
			if(gene== check):
				locus=info[0].split("..")
				left= int( locus[0] )
				right= int( locus[1] )
				window_ID= int(gene.split("_")[1])
				means=geneID_dic[gene]
				window_Dic[window_ID]= [left, right, means[0], means[1]]
				break
	file_ptt.close()
	
	combined_dic= combineWindowsDic(window_Dic)
	
	##This should only contain the CCNA_R0081 window information
	file_rnt= open(os.path.realpath("genome_files/CB15_SW/CP001340_2016_09_22_SW_25.rnt") )
	species = file_rnt.readline()
	rnt_Count= file_rnt.readline()
	header = file_rnt.readline()
	
	for line in file_rnt:
		for gene in geneIDs:
			if(line.find(gene) > -1):
				info=line.strip().split("\t")
				locus=info[0].split("..")
				left= int( locus[0] )
				right= int( locus[1] )
				means=geneID_dic[gene]
				combined_dic[gene]= [left, right, means[0], means[1], "", 1]
				break
	file_rnt.close()
	
	combine_keys= combined_dic.keys()
	
	annot_dic= make_gene_dic()
	
	new_dic={}
	## meanA[0], meanB[1], Log2Fold[2], gene_match[3]
	for cle in combine_keys:
		left= combined_dic[cle][0]
		right= combined_dic[cle][1]
		meanA= combined_dic[cle][2]
		meanB= combined_dic[cle][3]
		w_log= combined_dic[cle][4]
		log2= calculateLog2Fold(meanA, meanB)
		gene_match, strand= match_window_get_strand(annot_dic, left, right)
		gleft="-"
		gright="-"
		if(gene_match!="NA"):
			gene_info= annot_dic[gene_match]
			##left[0], right[1], strand[2], length[3], gene_name[4], product[5],
			gleft= gene_info[0]
			gright= gene_info[1] 
		new_dic[cle]= [left, right, meanA, meanB, log2, w_log, gene_match, strand, gleft, gright]
	
	return new_dic

def getWindowID_DeSeq(file_input_string):
	in_file=open(file_input_string)
	in_file.readline()
	out_window_dic={}
	for line in in_file:
		info=line.strip().split()
		window1 =info[0]
		baseA= float(info[2])
		baseB= float(info[3])
		out_window_dic[window1]= [baseA, baseB]
	return out_window_dic
	
def writeOutputDic(output_dic, file_name):
	file_out= open(file_name, "w")
	file_out.write("window_ID\tleft\tright\tPOS_mean\tNEG_mean\tLogFold2\tWindowLog\tgene_match\tstrand\tgene_left\tgene_right\n")
	
	out= output_dic.keys()
	for g in out:
		info= output_dic[g]
		info_string=""
		for i in range(len(info)-1):
			info_string+= str(info[i]) +"\t"
		info_string+= str(info[-1]) + "\n"
		file_out.write(g + "\t"+ info_string)
	file_out.close()
