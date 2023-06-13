
### Preclinical 2018 Imatinib Meta-Analysis Results ###

## Set up

library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import Data

df = read.csv("../../ReferenceFiles/PreclinicalMetaAnalysis/PreclinicalImatinib2018_20200805.csv",header=T,stringsAsFactors=F)

df$PMID = factor(df$PMID,levels=df$PMID[order(df$Concentration.nM)])

## Plot results

# Scatter plot
ggplot(df,aes(x=PMID,y=Concentration.nM/1e3))+theme_bw()+
  geom_point(size=3,color="#708090")+
  xlab("Preclinical Studies")+
  ylab("Lowest Efficacious Dose\nin in vitro System [uM]")+
  ggtitle("Imatinib Preclinical Meta-Analysis")+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  # scale_color_gradient(low="#66CC77",high="#FF4747",trans="log")+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=18,color="black",face="bold"),
    axis.text.x = element_blank(),
    legend.title = element_text(size=14,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=12,face="bold"),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  guides(color=F)

# ggsave("Imatinib2018Concentrations.png",width=8,height=6)

# Histogram (without outliers)

ggplot(df[df$Concentration.nM>=1e0&df$Concentration.nM<=1e6,],aes(x=Concentration.nM/1e3))+theme_minimal()+
  geom_histogram(bins=25,fill="#9fb7cf")+
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  ylab("Count")+
  theme(
    axis.title.y = element_text(size=24,color="black",face="bold"),
    # axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=20,color="black",face="bold",angle=90),
    axis.text.x = element_blank()
  )

# ggsave("Imatinib2018ConcentrationsHist.png",width=8,height=6)

