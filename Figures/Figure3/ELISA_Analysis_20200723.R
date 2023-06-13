
### Analyze ABL ELISA results ###

## Set up

library(reshape2)
library(ggplot2)
library(scales)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and parse ELISA results

raw.df = read.csv("../../ReferenceFiles/ELISAData_20191228.csv",header=T,stringsAsFactors=F)

# Standardize against vehicle control
stnd.df = cbind(raw.df[1:3],
                raw.df[4:12]/raw.df$vc)
                  

stnd.mlt = melt(stnd.df,id.vars = c("drug","serum","replicate"))
colnames(stnd.mlt)[4:5] = c("X.log.dose","absorbance")

stnd.mlt$log.dose = as.numeric(gsub("X","",stnd.mlt$X.log.dose))

## Aggregate data

agg.df = aggregate(absorbance~drug+serum+log.dose,stnd.mlt,function(x) c(mean = mean(x), sem = sd(x)/sqrt(2)))
agg.flt = do.call("data.frame",agg.df)

## Plot data

cmmn.thm = theme(
  plot.title = element_text(size=24,face="bold",hjust=0.5),
  axis.title = element_text(size=18,face="bold"),
  axis.text = element_text(size=14,color="black"),
  legend.title = element_text(size=14,face="bold"),
  legend.title.align = 0.5,
  legend.text = element_text(size=12,face="bold"),
  legend.position = c(.275,.2),
  legend.background = element_rect(color="black",fill="gray90"),
  legend.key = element_rect(color="black")
)

# Imatinib
ggplot(agg.flt[agg.flt$drug=="Imatinib",],aes(x=10^log.dose,y=absorbance.mean,color=serum))+theme_bw()+
  # geom_point(size=2.5)+
  geom_line(size=2.5)+
  xlab("Imatinib Dose [nM]")+ylab("Absorbance\nRelative to Untreated")+
  ggtitle("p-ABL ELISA: Imatinib")+
  geom_errorbar(aes(ymin=absorbance.mean-absorbance.sem,ymax=absorbance.mean+absorbance.sem),width=0.2,size=1.2)+
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_color_manual("Condition",labels=c("Without Serum Proteins","With Serum Proteins"),values=c("#999999","#785DAC"))+
  cmmn.thm

# ggsave("ELISA_Imatinib.png",width=6,height=4)


# Nilotinib
ggplot(agg.flt[agg.flt$drug=="Nilotinib",],aes(x=10^log.dose,y=absorbance.mean,color=serum))+theme_bw()+
  # geom_point(size=2.5)+
  geom_line(size=2.5)+
  xlab("Nilotinib Dose [nM]")+ylab("Absorbance\nRelative to Untreated")+
  ggtitle("p-ABL ELISA: Nilotinib")+
  geom_errorbar(aes(ymin=absorbance.mean-absorbance.sem,ymax=absorbance.mean+absorbance.sem),width=0.2,size=1.2)+
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_color_manual("Condition",labels=c("Without Serum Proteins","With Serum Proteins"),values=c("#999999","#785DAC"))+
  cmmn.thm

# ggsave("ELISA_Nilotinib.png",width=6,height=4)

