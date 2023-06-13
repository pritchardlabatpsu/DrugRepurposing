
### Plot Results of Imatinib Preclinical Meta-Analysis: Address PK ###

## Set up

library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and tally data

# Import
preclin.df = read.csv("../../ReferenceFiles/PreclinicalMetaAnalysis/PreclinicalImatinib2018_20200805.csv",header=T)

# Tally
tally.df = as.data.frame(table(preclin.df$PK.Addressed))
tally.df$Var1 = factor(tally.df$Var1,c("No","Uncorrected","Corrected"))

## Plot results

ggplot(tally.df,aes(x="",y=Freq,fill=Var1))+theme_void()+
  geom_bar(stat="identity",width=1,color="white")+
  coord_polar("y",start=0)+
  scale_fill_manual("",labels=c("Pharmacokinetics\nNot Addressed",
                                "Uncorrected Plasma\nConcentration Considered",
                                "Corrected Plasma\nConcentration Considered"),
                    values=c("#A0E3DC","#58979D","#63264D"))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    legend.title = element_text(size=14,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=12,face="bold")
    # legend.position = "bottom",
    # legend.spacing.x = unit(.5,'cm')
  )+
  guides(fill=F)

# ggsave("ImatinibPieChart.png",width=5.5,height=5.5)
