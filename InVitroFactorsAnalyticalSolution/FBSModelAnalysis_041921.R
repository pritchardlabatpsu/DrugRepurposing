
### Analyze results of biochemical kinetics simulations for off-target binding
### Scott Leighow - 04/19/21

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(reshape2)
library(ggplot2)
library(scales)
library(viridis)

## Import simulation results
vivo.df = read.csv("IC50vivo_041921.csv",header=F)
vitro.df = read.csv("IC50vitro_041921.csv",header=F)
vitroHSA.df = read.csv("IC50vitroHSA_041921.csv",header=F)

KB.fold = as.numeric(read.csv("K_BFoldIncreaseOverK_T.txt",sep=",",header=F))
KH.fold = as.numeric(read.csv("K_HFoldIncreaseOverK_T.txt",sep=",",header=F))

## Tidy data frames

colnames(vivo.df) = KH.fold
colnames(vitro.df) = KH.fold
colnames(vitroHSA.df) = KH.fold

vivo.df$KB.fold = KB.fold
vitro.df$KB.fold = KB.fold
vitroHSA.df$KB.fold = KB.fold

vivo.mlt = melt(vivo.df,id.vars = "KB.fold",variable.name = "KH.fold",value.name = "IC50")
vivo.mlt$KH.fold = as.numeric(as.character(vivo.mlt$KH.fold))
vitro.mlt = melt(vitro.df,id.vars = "KB.fold",variable.name = "KH.fold",value.name = "IC50")
vitro.mlt$KH.fold = as.numeric(as.character(vitro.mlt$KH.fold))
vitroHSA.mlt = melt(vitroHSA.df,id.vars = "KB.fold",variable.name = "KH.fold",value.name = "IC50")
vitroHSA.mlt$KH.fold = as.numeric(as.character(vitroHSA.mlt$KH.fold))

## Calculate fold-increase of IC50 over in vivo value for +HSA and -HSA conditions

fold.diff = data.frame(
  KB.fold = vivo.mlt$KB.fold,
  KH.fold = vivo.mlt$KH.fold,
  HSA.fold = vitroHSA.mlt$IC50/vivo.mlt$IC50,
  noHSA.fold = vitro.mlt$IC50/vivo.mlt$IC50
)

## Plot global results

ggplot(fold.diff,aes(x=KH.fold*0.01,y=KB.fold*0.01,fill=log10(HSA.fold)))+theme_bw()+
  geom_raster()+
  xlab("")+
  ylab("")+
  scale_x_log10(labels=trans_format('log10',math_format(10^.x)),breaks=c(1e-2,1e0,1e2))+
  scale_y_log10(labels=trans_format('log10',math_format(10^.x)),breaks=c(1e-2,1e0,1e2))+
  scale_fill_gradientn("EC50 Ratio",colors=RColorBrewer::brewer.pal(11,"Spectral"),limits=c(-4.5,3.5),
                       values=c(0,seq(0.4,0.55,length=4),seq(0.58,0.9,length=5),1))+
  theme(
    axis.text = element_text(size = 16,color="black"),
    legend.title = element_text(size=16,face="bold")
  )
# ggsave("IC50FoldDiff_HSA.png",width=6,height=4)

ggplot(fold.diff,aes(x=KH.fold*0.01,y=KB.fold*0.01,fill=log10(noHSA.fold)))+theme_bw()+
  geom_raster()+
  xlab("")+
  ylab("")+
  scale_x_log10(labels=trans_format('log10',math_format(10^.x)),breaks=c(1e-2,1e0,1e2))+
  scale_y_log10(labels=trans_format('log10',math_format(10^.x)),breaks=c(1e-2,1e0,1e2))+
  scale_fill_gradientn("EC50 Ratio",colors=RColorBrewer::brewer.pal(11,"Spectral"),limits=c(-4.5,3.5),
                       values=c(0,seq(0.4,0.55,length=4),seq(0.58,0.9,length=5),1))+
  theme(
    axis.text = element_text(size = 16,color="black"),
    legend.title = element_text(size=16,face="bold")
  )
# ggsave("IC50FoldDiff_noHSA.png",width=6,height=4)

## Import data for specific example (KB and KH 100x KT) and plot

ex.df = read.csv("ExampleCurve_041921.csv",header=F)
colnames(ex.df) = c("drug.conc","vivo","vitro","vitroHSA")

ex.mlt = melt(ex.df,id.vars="drug.conc",variable.name="condition",value.name="frc.bnd")

ggplot(ex.mlt[ex.mlt$drug.conc>=1e-3&ex.mlt$drug.conc<=1e1,],aes(x=drug.conc,y=frc.bnd,color=condition))+theme_light()+
  geom_line(size=2,alpha=0.8)+
  xlab("Total Drug Concentration")+
  ylab("Fraction of Target Bound")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_x_log10(labels=trans_format('log10',math_format(10^.x)))+
  ggtitle("Sample Dose Response Curve")+
  theme(
    plot.title = element_text(size=22,color="black",face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=16,color="black"),
  )
# ggsave("ExampleCurve.png",width=5.3,height=4)

  