
### Plotting Clinical Trial Summary ###

## Set up

library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import data

erl.df = read.csv("../../ReferenceFiles/RepurposingClinicalTrials/ErlotinibRepurposingTrialsSummary_20210214.csv",header=T)

## Plot data

ggplot(erl.df,aes(x=Approved,y=N,fill=Status))+theme_bw()+
  geom_bar(stat="identity")+
  xlab("")+ylab("Number of Unique Diseases")+
  ggtitle("Erlotinib Repurposing:\nClinical Trial Outcomes")+
  scale_x_discrete(labels=c("Not\nApproved","Approved"))+
  scale_fill_manual(values=c("#F53D00","#FDC835","#33673B"))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=18,color="black",face="bold"),
    axis.text.x = element_text(vjust=0.5),
    legend.title = element_text(size=14,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=12,face="bold")
  )

# ggsave("ErlotinibClinicalTrialResults.png",width=5.5,height=5.5)

## Alternative: Plot as pie chart

apprvd.df = aggregate(erl.df$N, by=list(Approved=erl.df$Approved), FUN=sum)
ggplot(apprvd.df, aes(x="", y=x, fill=Approved))+theme_void()+
  geom_bar(stat="identity", width=1,color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#D13428","#63C16C"))+
  guides(fill=F)
# ggsave("ErlotinibAnalysisPieChartOuter.png",width=8,height=8)

status.df = aggregate(erl.df$N, by=list(Status=erl.df$Status), FUN=sum)
ggplot(status.df, aes(x="", y=x, fill=Status))+theme_void()+
  geom_bar(stat="identity", width=1,color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#FF5B57","#FCD26B","#A6F4AC"))+
  guides(fill=F)
# ggsave("ErlotinibAnalysisPieChartInner.png",width=8,height=8)

ggplot(status.df, aes(x="", y=x))+theme_void()+
  geom_bar(stat="identity", width=1,fill=brewer.pal(10,"Spectral")[6])+
  coord_polar("y", start=0)
# ggsave("ErlotinibAnalysisPieChartRing.png",width=8,height=8)
