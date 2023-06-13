
### BCR-ABL IC50 Curves ###

## Set up

library(reshape2)
library(ggplot2)
library(scales)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load and clean data

IC50.df = read.csv("../../../ReferenceFiles/BCRABL_IC50s.csv",header=T,stringsAsFactors=F)
colnames(IC50.df)[1] = "genotype"

mutant = "T315I"

WT.df = IC50.df[IC50.df$genotype=="WT",]
mut.df = IC50.df[IC50.df$genotype==mutant,]

# Put data in tidy format
WT.df = WT.df[,grepl("n",colnames(WT.df))]
WT.mlt = melt(WT.df,id.vars="genotype")
colnames(WT.mlt) = c("genotype","dose.uM","frac.viable")
WT.mlt$dose.uM = as.numeric(gsub("n","",WT.mlt$dose.uM))

mut.df = mut.df[,grepl("n",colnames(mut.df))]
mut.mlt = melt(mut.df,id.vars="genotype")
colnames(mut.mlt) = c("genotype","dose.uM","frac.viable")
mut.mlt$dose.uM = as.numeric(gsub("n","",mut.mlt$dose.uM))

# Dosing parameters
Cave.free = 203/1e3
Cave.eff = 444/1e3

## QC Data

## Fit IC50 Curve to Data

# lm.fit = lm(log10(x)~y)
# a.start = as.numeric(10^predict(lm.fit,data.frame(y=0.5)))

WT.fit = nls(frac.viable~1/(1+(IC50/dose.uM)^hill),
             data = WT.mlt,
             start = list(IC50=.2,hill=2))
x.WT = 10^seq(-2,2,length.out=100)
y.WT = predict(WT.fit,newdata=list(dose.uM=x.WT))

mut.fit = nls(frac.viable~1/(1+(IC50/dose.uM)^hill),
             data = mut.mlt,
             start = list(IC50=.8,hill=2))
x.mut = 10^seq(-2,2,length.out=100)
y.mut = predict(mut.fit,newdata=list(dose.uM=x.mut))

## Plot Results

plot.df = rbind(
  data.frame(
    genotype = "WT",
    dose.uM = x.WT,
    frac.viab = y.WT
  ),
  data.frame(
    genotype = mutant,
    dose.uM = x.mut,
    frac.viab = y.mut
  )
)

ggplot(plot.df,aes(x=dose.uM,y=frac.viab,color=genotype))+theme_bw()+
  geom_vline(aes(xintercept=Cave.free),color="gray80",linetype="dashed",size=1.5)+
  geom_hline(aes(yintercept=0.5),color="gray80",linetype="dashed",size=1.5)+
  geom_line(size=3)+
  xlab("Dose [uM]")+
  ylab("Fraction Viable")+
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(0,1))+
  scale_color_manual("Genotype",values=c("#5AD891","#6f7a9b"),labels=c("WT","Mutant"))+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black"),
    axis.text.x = element_text(angle=45,vjust=0.5),
    legend.title = element_text(size=12,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=8,face="bold"),
    legend.position = c(.8,.8),
    legend.background = element_rect(color="black",fill="gray90"),
    legend.key = element_rect(color="black")
  )

# ggsave("BCRABL_IC50Curve.png",height=4,width=6)

# 50% viable line

ggplot(plot.df,aes(x=dose.uM,y=frac.viab,color=genotype))+theme_bw()+
  # geom_vline(aes(xintercept=Cave.free),color="gray80",linetype="dashed",size=1.5)+
  geom_hline(aes(yintercept=0.5),color="gray80",linetype="dashed",size=1.5)+
  geom_line(size=3)+
  xlab("Dose [uM]")+
  ylab("Fraction Viable")+
  scale_x_continuous(trans=log10_trans(),
                     # breaks=trans_breaks("log10", function(x) 10^x),
                     breaks = c(1e-2,1e0,1e2),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(0,1),labels=c("0.0","1.0"))+
  scale_color_manual("Genotype",values=c("#5AD891","#6f7a9b"),labels=c("WT","Mutant"))+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black"),
    axis.text.x = element_text(angle=45,vjust=0.5),
    legend.title = element_text(size=12,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=8,face="bold"),
    legend.position = c(.8,.8),
    legend.background = element_rect(color="black",fill="gray90"),
    legend.key = element_rect(color="black"),
    panel.grid = element_blank()
  )

# ggsave("BCRABL_IC50CurveHline.png",height=4,width=6)

# free Cave line

ggplot(plot.df,aes(x=dose.uM,y=frac.viab,color=genotype))+theme_bw()+
  geom_vline(aes(xintercept=Cave.free),color="gray80",linetype="dashed",size=1.5)+
  # geom_hline(aes(yintercept=0.5),color="gray80",linetype="dashed",size=1.5)+
  geom_line(size=3)+
  xlab("Dose [uM]")+
  ylab("Fraction Viable")+
  scale_x_continuous(trans=log10_trans(),
                     # breaks=trans_breaks("log10", function(x) 10^x),
                     breaks = c(1e-2,1e0,1e2),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(0,1),labels=c("0.0","1.0"))+
  scale_color_manual("Genotype",values=c("#5AD891","#6f7a9b"),labels=c("WT","Mutant"))+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black"),
    axis.text.x = element_text(angle=45,vjust=0.5),
    legend.title = element_text(size=12,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=8,face="bold"),
    legend.position = c(.8,.8),
    legend.background = element_rect(color="black",fill="gray90"),
    legend.key = element_rect(color="black"),
    panel.grid = element_blank()
  )

# ggsave("BCRABL_IC50CurveVline.png",height=4,width=6)

## Plot serum shift example

y.WT.shft = predict(WT.fit,newdata=list(dose.uM=x.WT/10))

plot2.df = rbind(
  data.frame(
    condition = "Without Serum Proteins",
    dose.uM = x.WT,
    frac.viab = y.WT
  ),
  data.frame(
    condition = "With Serum Proteins",
    dose.uM = x.mut,
    frac.viab = y.WT.shft
  )
)

ggplot(plot2.df,aes(x=dose.uM,y=frac.viab,color=condition))+theme_bw()+
  geom_hline(aes(yintercept=0.5),color="gray80",linetype="dashed",size=1.5)+
  # geom_vline(aes(xintercept=Cave.eff),color="gray80",linetype="dashed",size=1.5)+
  # geom_vline(aes(xintercept=Cave.eff/10),color="gray80",linetype="dashed",size=1.5)+
  geom_line(size=3)+
  xlab("Dose [uM]")+
  ylab("Fraction Viable")+
  scale_x_continuous(trans=log10_trans(),
                     breaks = c(1e-2,1e2),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(0,1))+
  scale_color_manual("Condition",values=c("#785DAC","#999999"))+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black"),
    axis.text.x = element_text(angle=45,vjust=0.5),
    legend.title = element_text(size=16,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=12,face="bold"),
    legend.position = c(.75,.8),
    legend.background = element_rect(color="black",fill="gray90"),
    legend.key = element_rect(color="black")
  )

# ggsave("SerumShiftIC50a.png",height=4,width=6)

ggplot(plot2.df[plot2.df$condition=="Without Serum Proteins",],aes(x=dose.uM,y=frac.viab,color=condition))+theme_bw()+
  # geom_hline(aes(yintercept=0.5),color="gray80",linetype="dashed",size=1.5)+
  geom_vline(aes(xintercept=Cave.eff),color="#79CCCC",linetype="dashed",size=1.5)+
  geom_vline(aes(xintercept=Cave.eff/10),color="gray80",linetype="dashed",size=1.5)+
  geom_line(size=3)+
  xlab("Dose [uM]")+
  ylab("Fraction Viable")+
  scale_x_continuous(trans=log10_trans(),
                     breaks = c(1e-2,1e2),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(0,1))+
  scale_color_manual("Condition",values=c("#999999"))+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black"),
    axis.text.x = element_text(angle=45,vjust=0.5),
    legend.title = element_text(size=16,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=12,face="bold"),
    legend.position = c(.75,.8),
    legend.background = element_rect(color="black",fill="gray90"),
    legend.key = element_rect(color="black")
  )

# ggsave("SerumShiftIC50b.png",height=4,width=6)

ggplot(plot2.df[plot2.df$condition=="Without Serum Proteins",],aes(x=dose.uM,y=frac.viab,color=condition))+theme_bw()+
  # geom_hline(aes(yintercept=0.5),color="gray80",linetype="dashed",size=1.5)+
  geom_vline(aes(xintercept=Cave.eff/3*10),color="gray80",linetype="dashed",size=1.5)+
  geom_vline(aes(xintercept=Cave.eff/3),color="#79CCCC",linetype="dashed",size=1.5)+
  geom_line(size=3)+
  xlab("Dose [uM]")+
  ylab("Fraction Viable\n")+
  scale_x_continuous(trans=log10_trans(),
                     breaks = c(1e-2,1e2),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(0,1))+
  scale_color_manual("Condition",values=c("#999999"))+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black"),
    axis.text.x = element_text(angle=45,vjust=0.5),
    legend.title = element_text(size=16,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=12,face="bold"),
    legend.position = c(.75,.8),
    legend.background = element_rect(color="black",fill="gray90"),
    legend.key = element_rect(color="black")
  )

# ggsave("SerumShiftIC50c.png",height=4,width=6)

