
### Heat Maps of Classifier Variables ###

## Set up

library(ggplot2)
library(RColorBrewer)
library(scales)

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


# Free alpha
ggplot(alpha.df,aes(x=drug,y=genotype,fill=free.alpha))+theme_bw()+
  geom_raster()+
  xlab("Drugs")+ylab("Genotype")+
  ggtitle("Drug Kill Rate at Free Cave")+
  scale_fill_gradientn("Free\nAlpha",colors=brewer.pal(10,"RdYlBu"))+
  cmmn.thm

# ggsave("FreeAlpha_HeatMap_20200722.png",width=6,height=6)

# IC50 fold change
ggplot(alpha.df,aes(x=drug,y=genotype,fill=IC50.fld.chng))+theme_bw()+
  geom_raster()+
  xlab("Drugs")+ylab("Genotype")+
  ggtitle("Fold Change in IC50 over WT")+
  scale_fill_gradientn("IC50\nFold-Change",colors=rev(brewer.pal(10,"RdYlBu")),
                       trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x)))+
  cmmn.thm

# ggsave("IC50FoldChange_HeatMap_20200722.png",width=6,height=6)

