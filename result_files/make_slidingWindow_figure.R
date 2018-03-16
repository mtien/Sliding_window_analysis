setwd("/Users/mtien/Documents/Crosson/Manuscripts/Elife/Figure_source_data/Figure5_sourcedata2/")
library(ggplot2)
library(ggrepel)
library(reshape)
library(extrafont)
library(cowplot)

temp<-read.delim("POSvsNEG_105.txt")
temp<-temp[temp$pval<.1, ]
temp<-temp[log10(temp$baseMeanA)>2.5, ]
temp<-data.frame(id=temp$id, Xval=log(-log10(temp$pval)), Yval= log10(temp$baseMeanA), Sig="InSig")

temp2<-read.delim("POSvsNEG_p10_BMA1000_POS_windowInfo_RhopperCongruent.txt")
temp2<-data.frame(id=temp2$id, Xval=log(-log10(temp2$pval)), Yval= log10(temp2$baseMeanA), Sig="Congruent")

temp3<-read.delim("POSvsNEG_p10_BMA1000_POS_windowInfo_RhopperNotCongruent.txt")
temp3<-data.frame(id=temp3$id, Xval=log(-log10(temp3$pval)), Yval= log10(temp3$baseMeanA), Sig="Not Congruent")

temp_sig<-rbind(temp2, temp3)
temp<- temp[!(temp$id %in% temp_sig$id),]
data_frame<- rbind(temp_sig,temp)

geom.text.size = 3*.7
theme.size = (14/5) * geom.text.size
p<-ggplot(data_frame, aes(x=Xval, y=Yval, color=Sig)) + 
  geom_point() +
  geom_hline(yintercept=3, color="red", lty=2) + 
  geom_vline(xintercept = log(-log10(0.1)), color="red", lty=2) +
  lims(y=c(2.5,6) )+ guides(color=FALSE) +
  scale_color_manual(values=c("blue","orange","black")) +
  labs(x="log(-log10(qValue))", y="gsrN(37)PP7hp Expression")+
  theme(axis.text = element_text(size = theme.size), axis.title = element_text(size = theme.size)) 

ggsave("Figure4_SlidingWindow_ggplot.pdf", p, width=80, height=56.5, units = c("mm"), dpi = 300)

