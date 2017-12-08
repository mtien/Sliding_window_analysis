setwd("~/Desktop/BigDisk/Sliding_window/Github")
library("DESeq")
countTable= read.delim( "data_files/RPKM_compiled_slidingWindow_removeHighVariantWindows.txt", header=TRUE, row.names=1 )
conditionTable = data.frame(row.names = colnames(countTable), condition = c("total","total","total","total", "POS","POS","POS", "NEG", "NEG") )

cds= newCountDataSet( countTable, conditionTable$condition )
cds = estimateSizeFactors( cds )
cds= estimateDispersions( cds )

res_No_WvR = nbinomTest( cds, "POS", "NEG" )
resSig_No_WvR = res_No_WvR[ res_No_WvR$padj < 1.05, ]
write.table( resSig_No_WvR[ order(resSig_No_WvR$pval), ], file= "data_files/POSvsNEG_105.txt", quote= FALSE, sep="\t", row.names=FALSE )
