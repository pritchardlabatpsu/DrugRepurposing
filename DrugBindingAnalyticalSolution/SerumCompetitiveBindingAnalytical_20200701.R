
### Plot analytical solution to competitive binding equations ###

## Set up

library(reshape2)
library(scales)
library(ggplot2)
library(viridis)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import results
htmp.df = read.csv("DrugNotBoundToSerum_20200701.csv",header=F)
colnames(htmp.df) = c("Kd.T","Kd.S",       # parameters
                      "DnS.wc","DnS.woc",  # fraction not bound to serum (wc=with competition; woc=without)
                      "DnS.ratio")         # ratio of DnS.wc to DnS.woc

## Plot results

ggplot(htmp.df,aes(x=as.factor(log10(Kd.S)),y=as.factor(log10(Kd.T)),fill=log10(DnS.ratio)))+theme_minimal()+
  geom_tile()+
  xlab("Serum Binding Affinity [uM]")+
  ylab("Target Binding Affinity [uM]")+
  ggtitle("Drug Not Bound to Serum: Ratio of\nCompetitive to Noncompetitive Binding")+
  scale_x_discrete(breaks=c(-3,2),labels=expression(10^-3,10^2))+
  scale_y_discrete(breaks=c(-5,1),labels=expression(10^-5,10^1))+
  scale_fill_gradientn("Ratio",
                       colors=c(viridis(10),rep(viridis(10)[10],3)),
                       # values = c(0,0.05,0.1,2,max(log10(htmp.df$DnS.ratio))),
                       limits=c(0,max(log10(htmp.df$DnS.ratio))),
                       breaks=c(0,2,4),
                       labels=expression(10^0,10^2,10^4))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=14,color="black"),
    legend.title = element_text(size=12,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=14,face="bold")
  )

# ggsave("SerumCompetitiveBindingHeatMap_20200702.png",width=8*.75,height=6*.75)
