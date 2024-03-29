
### Heat Maps of Classifier Variables ###

## Set up

library(ggplot2)
library(RColorBrewer)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import data

alpha.df = read.csv("../../../ReferenceFiles/AlphaValues_True_positive.csv",header=T,stringsAsFactors=F)

## Intermediate Calculations

# Calculate IC50 fold change
drugs = unique(alpha.df$drug)

for (drug in drugs) {
  WT.IC50 = alpha.df$IC50[alpha.df$genotype=="WT"&alpha.df$drug==drug]
  alpha.df$IC50.fld.chng[alpha.df$drug==drug] = alpha.df$IC50[alpha.df$drug==drug]/WT.IC50
}

## Plot heat maps for classifiers

# Reorder drugs
alpha.df$drug = factor(alpha.df$drug,levels=c("Imatinib","Nilotinib","Dasatinib","Bosutinib","Ponatinib"))

# Reorder genotypes
vars = unique(alpha.df$genotype)
vars = vars[vars!="WT"]
vars = vars[order(substr(vars,2,5))]
alpha.df$genotype = factor(alpha.df$genotype,levels=rev(c("WT",vars)))

# Common theme
cmmn.thm = theme(
  plot.title = element_text(size=20,face="bold",hjust=0.5),
  axis.title = element_text(size=18,face="bold"),
  axis.text = element_text(size=12,color="black",face="bold"),
  axis.text.x = element_text(angle=15,vjust=0.5),
  legend.title = element_text(size=14,face="bold"),
  legend.title.align = 0.5,
  legend.text = element_text(size=14,face="bold")
)


# Effective alpha
ggplot(alpha.df,aes(x=drug,y=genotype,fill=eff.alpha))+theme_bw()+
  geom_raster()+
  xlab("Drugs")+ylab("Genotype")+
  ggtitle("Drug Kill Rate at Effective Cave")+
  scale_fill_gradientn("Effective\nAlpha",colors=brewer.pal(10,"RdYlBu"))+
  cmmn.thm

# ggsave("EffectiveAlpha_HeatMap_20200722.png",width=6,height=6)
