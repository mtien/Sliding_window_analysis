This is deposit of several python scripts, written by Matthew Z. Tien,

Reference:
Matthew Z. Tien, Aretha Fiebig, Sean Crosson (2017).
Gene network analysis identifies a central post-transcriptional regulator of cellular stress survival
bioRxiv 212902; doi: https://doi.org/10.1101/212902 

The most up-to-date version of this software is
available at https://github.com/mtien/Sliding_window_analysis.

Rockhopper analysis:

	Rockhopper software package generates a "transcripts.txt" tab-delimited file when the "verbose output"
	option is turned on. The script, parse_Rhopper_transcript_file.py, takes the criterion outlined 
	in the Materials and Methods section from bioRxiv 212902; doi: https://doi.org/10.1101/212902
	and creates a "transcripts_criterion.txt" file. This output file is important when comparing the 
	Sliding window analysis with a standard RNA-Seq analysis approach.
	
	The R script, make_qval_vs_abundance_Rhopper.R, can use these two files to generate figure 5C in
	bioRxiv 212902; doi: https://doi.org/10.1101/212902 

Sliding Window analysis:

	Several libraries were constructed to analyze the PP7-purification total RNA-Seq read data.
	These libraries include: assign_window.py, check_bowtie_alignment.py, and synthesize_information.py .

		assign_window.py contains a series of methods that will take window_IDs 
		from the genome_files directory and map them to proximal genes around the window_ID of
		interest. It also contains methods to combine adjacent sliding windows. 

		check_bowtie_alignment.py contains a series of methods to parse a bowtie alignment file (".hits") 
		that has already been parsed to contain the following FLAGS: 0, 256, 16, 272. Unix command line used to generate a bowtie "hits" file:
		awk '{split($0,arr,\"\\t\"); if(arr[2]==\"0\" || arr[2]==\"256\" || arr[2]==\"16\" || arr[2]==\"272\") print $0}' MZT_sense1.alignments > MZT_sense1.alignments.hits
		where MZT_sense1.alignments correspond to the bowtie output file generated from running the EDGE-pro software package.

		synthesize_information.py is a variant of the assign_window library, but helps incorporate the
		information from the check_bowtie_alignment library.
	
	A series of python scripts utilized the libraries to perform the sliding window analysis.

		As described in the Materials and Methods section of RNA-seq analysis of mRNAs that co-elute with GsrN doi: 
		https://doi.org/10.1101/212902 RNA-seq analysis of mRNAs that co-elute with GsrN, removal occurred before analysis
		by DESeq in order to decrease the False Positive Rate and to balance the read density between the PP7 purifications.
		The script that corresponds to this process is called remove_high_variant_windows.py. This script takes in the 
		"RPKM_compiled_slidingWindow.txt" file. This file is zipped on the GitHub under the data_files folder.
	
		After removal of inconsistent windows, a DESeq script is run called, deSeq.R. This generates a table of all sliding windows
		and their significance as judged by the DESeq software package.

		After the DESeq estimates the significance of each sliding window, the assign analysis.py will take the DESeq output file and do
		a similar analysis to the Rockhopper analysis script, parse_Rhopper_transcript_file.py.
		assign_analysis utilizes the three libraries: assign_window, check_bowtie_alignment, and
		synthesize_information to create the results file ("Sliding_Window_analysis_results.txt") of the sliding window analysis.
		The script will generate several intermediate files and produces the resulting file.

Combined analysis:

	Once the Rockhopper analysis and sliding window analysis have generated their final result files, several scripts can
	be used to compare the results of each analysis. The first script that should be run is get_window_information.py, 
	which will take one of the intermediate files from the sliding window analysis and breaks the windows down into a
	new "windowInfo.txt" file. 
	
	The script corroborate_SW_Rhopper.py takes the "windowInfo.txt" file and the "transcripts_criterion.txt" file and generates
	several files that show which genes overlap in these two separate analyses. The files "RhopperCongruent.txt" and 
	"RhopperNotCongruent.txt" were used to generate Figure 5D in doi: https://doi.org/10.1101/212902 using the R script,
	make_sliding_window_figure.R.
	
	The final script, make_final_results_table.py, takes the two output files of the sliding window analysis and 
	Rockhopper analysis and combines them into a "Final_result_table.txt", Table 1 in
	doi: https://doi.org/10.1101/212902.

IntaRNA analysis:

	FASTA_methods is a light weight python library to retrieve sequences from Caulobacter's
	genome based on the input of gene coordinates.
	
	The script, get_FASTA_sequences.py, utilizes the FASTA_methods library to retrieve the sequences identified
	in the two output files of the sliding window analysis and Rockhopper analysis. The final FASTA output file
	can then be run with the IntaRNA software suite.
	
	The final script parse_intaRNA.py takes the csv file from the IntaRNA online package and creates several files
	highlighting where GsrN most likely associates with the inputted FASTA file and what part of GsrN is most used
	in interacting with it's binding partners.
