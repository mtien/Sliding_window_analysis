#!usr/bin/python

import assign_window

def match_window(g_dic, left, right, strand):
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
			position, magnitude= assign_window.determine_position(g_strand, left_g, right_g, left, right)
			pos_dic[p]=[position, magnitude, g_strand, left_g, right_g]
		
		gene_internal=""
		internal_count= 0
		pd= pos_dic.keys()
		sense="S"
		for p in pd:
			if(pos_dic[p][0] == "I"):
				gene_internal= p
				internal_count+=1
				if(pos_dic[gene_internal][2] != strand):
					sense="A"
		
		if(internal_count == 1):
			return gene_internal, "I", sense
		elif(internal_count < 1):
			min_mag=250
			for p in pd:
				mag= pos_dic[p][1]
				if(mag < min_mag):
					min_mag = mag
					gene_internal=p
					if(pos_dic[gene_internal][2] != strand):
						sense="A"
						
			return gene_internal, pos_dic[gene_internal][0], sense
		else:
			gene_internal=""
			for p in pd:
				if(pos_dic[p][0] == "I"):
					gene_internal+=p +","
			return gene_internal, "O", sense
			
	elif( len(proximity) == 1 ):
		gene_ID= proximity[0]
		##left, right, strand, length, gene, product
		entries= g_dic[gene_ID]
		g_strand= entries[2]
		left_g= int(entries[0])
		right_g= int(entries[1])
		position="" 
		magnitude=""
		sense="S"
		position, magnitude= assign_window.determine_position(g_strand, left_g, right_g, left, right)
		if(g_strand != strand):
			sense="A"
		return gene_ID, position, sense 
	else:
		## NO GENES in proximity
		return "NA", "NA", "A"
		
def synthesize_information(file_in_name, file_out_name):
	gene_dic=assign_window.make_gene_dic()
	
	##window[0], left[1], right[2], posmean[3], negmean[4]. log2Fold[5],
	##windowLog[6], geneMatch[7], geneStrand[8], gene_left[9], gene_right[10]
	##btieStrand[15]
	file_in= open(file_in_name)
	file_in.readline()
	
	##remove rRNA matches
	gene_not=["CCNA_R0069", "CCNA_R0066", "CCNA_R0065", "CCNA_R0087", "CCNA_R0084", "CCNA_R0083"]
	
	##key combined window
	##left[0], right[1], btieStrand[2], PosMean[3], NegMean[4], length[5], 
	##associated gene(s)[6], categoryA[8] (I= internal, U= upstream, O=mapped to operon, D= downstream), categoryB[9] A=antisense S=sense
	##annotated gene left[10], annotated gene right[11], strand of annotated_gene [12]
	window_dic={}
	for line in file_in:
		info= line.strip().split("\t")
		if(info[7] not in gene_not):
			window_ID=info[0]
			w_left= info[1]
			w_right= info[2]
			btie_strand= info[15]
			pos_mean=info[3]
			neg_mean=info[4]
			log2fold=info[5]
			windowLog=info[6]
			geneMatch, catA, catB= match_window(gene_dic, int(w_left), int(w_right), btie_strand)
			min_left=4000000000
			max_right=0
			g_strand=""
			if(catA=="O"):
				genes= geneMatch.split(",")[:-1]
				##left[0], right[1], strand[2], length[3], gene_name[4], product[5]
				for g in genes:
					if(int(min_left) > int(gene_dic[g][0])):
						min_left=gene_dic[g][0]
					if(int(max_right) < int(gene_dic[g][1])):
						max_right= gene_dic[g][1]
					if(g_strand==""):
						g_strand= gene_dic[g][2]
					elif(g_strand != gene_dic[g][2]):
						g_strand= "MISMATCH"
			elif(geneMatch!="NA"):
				min_left=gene_dic[geneMatch][0]
				max_right= gene_dic[geneMatch][1]
				g_strand= gene_dic[geneMatch][2]
			else:
				min_left=w_left
				max_right=w_right
				g_strand=btie_strand
			##Calculate distance between window start and gene start
			TSS= int(w_left)-int(min_left)
			term= int(w_right) - int(max_right)
			if( g_strand == "-"):
				TSS= int(max_right)-int(w_right)
			if( g_strand == "-"):
				term= int(min_left)- int(w_left)
			window_dic[window_ID]=[w_left, w_right, btie_strand, pos_mean, neg_mean, log2fold, windowLog, geneMatch, catA, catB, min_left, max_right, g_strand, TSS, term]
	file_in.close()
	
	file_out= open(file_out_name, "w")
	file_out.write("window_ID\tleft_window\tright_window\tbtie_strand\tPOS_mean\tNEG_mean\tlog2Fold\tWindow_log\tGeneMatches\tCategory\tSense\tmin_left\tmax_right\tgene_strand\tTSS_distance\tterm_distance\n")
	windows= window_dic.keys()
	
	for w in windows:
		file_out.write(w)
		entries=window_dic[w]
		for e in entries:
			file_out.write("\t" + str(e) )
		file_out.write("\n")
	file_out.close()
