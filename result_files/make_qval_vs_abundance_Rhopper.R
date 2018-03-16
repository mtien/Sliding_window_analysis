setwd("/Users/mtien/Documents/Crosson/Manuscripts/Elife/Figure_source_data/Figure5_sourcetable1/")

library(ggplot2)
library(ggrepel)
library(reshape)
library(extrafont)
library(cowplot)

temp<-read.delim("NC_999999_transcripts.txt")
temp_all<-temp[ log(-log10(temp$qValue.POS.vs.NEG))>0,]
temp_all$qValue.POS.vs.NEG[temp_all$qValue.POS.vs.NEG ==0] <- 1e-260

temp_A<- temp_all[ temp_all$Expression.POS > temp_all$Expression.NEG, ]
temp_A<- temp_A[ temp_A$Expression.POS>1000,]
temp_A<- temp_A[ temp_A$qValue.POS.vs.NEG < .05, ]
temp_A_synonyms<- as.character(temp_A$Synonym)
temp_A_names<- as.character(temp_A$Name)

A_x<-log(-log10(temp_A$qValue.POS.vs.NEG))
A_y<-log10(temp_A$Expression.POS)
A_labels<- c()
for(i in seq(length(temp_A_synonyms)))
{
  if("-" == b[i]){A_labels<-c(A_labels,substring(temp_A_synonyms[i],6))}
  else{A_labels<-c(A_labels,temp_A_names[i])}
}
A_df<- data.frame(Sig="Sig", Xval=A_x, Yval=A_y, Label=A_labels)

temp_B<- temp_all[!(temp_all$Synonym %in% temp_A$Synonym),]
B_x<-log(-log10(temp_B$qValue.POS.vs.NEG))
B_y<-log10(temp_B$Expression.POS)
temp_B$Synonym<- as.character(temp_B$Synonym)
B_labels<- temp_B$Synonym
B_df<- data.frame(Sig="Not", Xval=B_x, Yval=B_y, Label=B_labels)

data_frame<-rbind(A_df, B_df)

geom.text.size = 3*.7
theme.size = (14/5) * geom.text.size

p<- ggplot(data_frame, aes(x=Xval, y=Yval, color=Sig)) + 
  geom_point() +
  geom_hline(yintercept=3, color="red", lty=2) + 
  geom_vline(xintercept = log(-log10(0.05)), color="red", lty=2) +
  lims(y=c(2,5) )+ guides(color=FALSE) +
  scale_color_manual(values=c("blue","black")) +
  geom_text_repel(aes(label=as.character(Label)), data=A_df, color="blue", size=geom.text.size) +
  labs(x="log(-log10(qValue))", y="gsrN(37)PP7hp Expression")+
  theme(axis.text = element_text(size = theme.size), axis.title = element_text(size = theme.size)) 
  
ggsave("Figure4_Rockhopper_ggplot.pdf", p, width=80, height=56.5, units = c("mm"), dpi = 300)



