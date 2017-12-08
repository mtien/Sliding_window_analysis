#!usr/bin/python

import sys, os
import assign_window
import check_bowtie_alignment
import synthesize_information

deSeq_output_file_name= "data_files/POSvsNEG_105.txt"
file_in= open(deSeq_output_file_name)
header= file_in.readline()
applied_criterium_file_name= "data_files/POSvsNEG_p10_BMA1000_POS.txt"
file_out= open(applied_criterium_file_name, "w")
file_out.write(header)
for line in file_in:
	info=line.strip().split()
	if( (line.find("NA\tNA") < 0)  and (float(info[2]) >1000) and (float(info[6]) < .10) and (float(info[2]) > float(info[3])) ):
		file_out.write(line)

file_in.close()
file_out.close()


window_dic=  assign_window.getWindowID_DeSeq(applied_criterium_file_name)
output_dic= assign_window.getCombinedWindowInfo(window_dic)
assigned_combined_file_name= "data_files/POSvsNEG_p10_BMA1000_POS_assignedCombinedWindows.txt"
assign_window.writeOutputDic(output_dic, assigned_combined_file_name)


window_position37, window_positionsPP7 = check_bowtie_alignment.getWindowDic(assigned_combined_file_name)
MZT1_dic= check_bowtie_alignment.parseBtieFile("alignment_file/MZT_sense1.alignments.hits", window_position37)
btie_info_file_name= "data_files/POSvsNEG_p10_BMA1000_POS_assignedCombinedWindows_btieInfo.txt"
check_bowtie_alignment.writeBtieInformation(MZT1_dic, assigned_combined_file_name, btie_info_file_name)

compiled_information_name= "result_files/Sliding_Window_analysis_results.txt"
synthesize_information.synthesize_information(btie_info_file_name, compiled_information_name)
