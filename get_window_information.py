#!usr/bin/python

import sys, os
import assign_window

if __name__ == "__main__":
	dic=assign_window.make_gene_dic()
	
	file_in= open("data_files/POSvsNEG_p10_BMA1000_POS.txt")
	header= file_in.readline().strip()+ "\tproximal_genes\tleft\tright\n"
	
	file_out= open("data_files/POSvsNEG_p10_BMA1000_POS_windowInfo.txt", "w")
	file_out.write(header)
	for line in file_in:
		info=line.strip().split()
		window1 =info[0]
		
		window_info=assign_window.find_window(window1)
		proxy= assign_window.match_window(dic, window_info)
		proxy_string=""
		for p in proxy:
			proxy_string+=p +","
		file_out.write(line.strip()+ "\t" + proxy_string + "\t" + window_info[0] + "\t" + window_info[1] +"\n")

	file_out.close()
	file_in.close()
