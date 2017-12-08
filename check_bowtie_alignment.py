#!usr/bin/python

import sys, os

def find_window(window, file_window_name="genome_files/CB15_SW/CP001340_2016_09_22_SW_25.ptt"):
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

def getWindowLeftRight(window):
	gsrN37=[]
	PP7=[]
	##find window returns array of 3
	if(window!="R0081"):
		gsrN37= find_window(window, file_window_name="genome_files/G12_SW/CP001340_2016_09_22_G12_SW_25_sense.ptt")
		PP7= find_window(window, file_window_name="genome_files/G17_SW/CP001340_2016_09_22_G17_SW_25_sense.ptt")
	else:
		gsrN37= find_window(window, file_window_name="genome_files/G12_SW/CP001340_2016_09_22_G12_SW_25_sense.rnt")
		PP7= find_window(window, file_window_name="genome_files/G17_SW/CP001340_2016_09_22_G17_SW_25_sense.rnt")
	return gsrN37[:2] + PP7[:2]
	
def getWindowDic(file_name):
	file_in= open(file_name)
	header= file_in.readline()
	##window_ID[0], left[1], right[2], POS_mean[3], NEG_mean[4], LogFold2[5], WindowLog[6], gene_match[7], strand[8], gene_left[9], gene_right[10]
	w_dic37={}
	w_dicPP7={}
	for line in file_in:
		info=line.strip().split()
		window_ID=info[0]
		window_left=int(info[1])
		window_right=int(info[2])
		if(float(info[5]) > 1.0):
			windowInfo=[]
			if( window_right-window_left > 25 and window_ID!="CCNA_R0081"):
				window_log=info[6].split(",")
				last_window= window_log[-2]
				firstWindow= getWindowLeftRight(window_ID)
				lastWindow= getWindowLeftRight(last_window)
				windowInfo=[ firstWindow[0], lastWindow[1], firstWindow[2], lastWindow[3] ]
			else:
				windowInfo= getWindowLeftRight(window_ID)
			w_dic37[window_ID]= windowInfo[:2]
			w_dicPP7[window_ID]= windowInfo[2:]
	file_in.close()
	return w_dic37, w_dicPP7
	
def parseBtieFile(btie_path, window_flag_dic):
	file_in= open(btie_path)
	
	w_keys=window_flag_dic.keys()
	flag_dic={key: [0,0,0,0] for key in w_keys}
	print "start bowtie assignment"
	for line in file_in:
		info= line.strip().split("\t")
		pos= int(info[3])
		for w in w_keys:
			window_left= int(window_flag_dic[w][0])
			window_right= int(window_flag_dic[w][1])
			if(pos >= window_left and pos <=window_right):
				flag= info[1]
				if(flag=="0"):
					flag_dic[w][0]+=1
				elif(flag=="16"):
					flag_dic[w][2]+=1
				elif(flag=="256"):
					flag_dic[w][1]+=1
				else:
					flag_dic[w][3]+=1
	file_in.close()
	print "finish bowtie assignment"
	return flag_dic
	
def writeBtieInformation(btie_dic, file_in_name, file_out_name):
	file_in= open(file_in_name)
	header= file_in.readline().strip() + "0_flag\t256_flag\t16_flag\t272_flag\tBtieStrand\n"
	
	windows= btie_dic.keys()
	
	file_out= open(file_out_name, "w")
	file_out.write(header)
	for line in file_in:
		info=line.strip().split("\t")
		window_id= info[0]
		output_string="-\t-\t-\t-\t-"
		if(window_id in windows):
			strand= "+"
			part0= str(btie_dic[window_id][0])
			part16= str(btie_dic[window_id][2])
			part256= str(btie_dic[window_id][1])
			part272= str(btie_dic[window_id][3])
			if(int(part0)+int(part256) < (int(part16)+int(part272)) ):
				strand="-"
			output_string= part0 +"\t"+ part256 +"\t"+ part16+"\t"+ part272 + "\t" +strand	 
		file_out.write( line.strip() )
		file_out.write( "\t" + output_string + "\n")
	file_in.close()
	file_out.close()
	
