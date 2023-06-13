
### Plot cell viability for range of HSA/AAG concentrations and drug conditions ###

## Set up

library(rstudioapi)
library(reshape2)
library(ggplot2)
library(scales)

rm(list=ls())
setwd(dirname(getActiveDocumentContext()$path))

## Loop through csv files and plot

no.shft.df = read.csv("NoShiftViabilities_20200723.csv",header=T,stringsAsFactors=F)

csv.files = grep("\\.csv$",list.files("TitrationData/"),value=T)

for (i in 1:length(csv.files)) {
  
  # Import data
  csv.i = csv.files[i]
  df.i = read.csv(paste("TitrationData/",csv.i,sep=""),header=T)
  
  # Clean up data
  nrep = ifelse(grepl("No_Drug",csv.i),6,2) # define number of replicates for each condition
  
  colnames(df.i) = c("Protein",
                     paste("HSA.AAG",1:nrep,sep="_"),
                     paste("HSA.Only",1:nrep,sep="_"),
                     paste("AAG.Only",1:nrep,sep="_"))
  
  # Organize into tidy format, by replicate
  df.mlt = melt(df.i,id.vars="Protein")
  
  df.mlt$Condition = unlist(lapply(strsplit(as.character(df.mlt$variable),"_"),'[[',1))
  df.mlt$Replicate = unlist(lapply(strsplit(as.character(df.mlt$variable),"_"),'[[',2))
  
  # Aggregate data
  df.agg = aggregate(df.mlt$value,
                     list(Protein=df.mlt$Protein,Condition=df.mlt$Condition),
                     FUN = function(x) c(mean=mean(x),sd=sd(x)))
  df.agg2 = do.call("data.frame",df.agg)
  
  # Get baseline (no serum protein) viability
  no.shft = no.shft.df$no.shift.viability[no.shft.df$condition==gsub("\\.csv","",csv.i)]
  df.agg2 = rbind(
    df.agg2,
    data.frame(
      Protein = c(min(df.agg2$Protein),max(df.agg2$Protein)),
      Condition = c("None","None"),
      x.mean = c(no.shft,no.shft),
      x.sd = c(NA,NA)
    )
  )
  
  # Plot data
  plot.i = gsub("\\.csv","",csv.i)
  title.i = gsub("_"," ",plot.i)
  title.i = gsub("0nM","0 nM",title.i)
  
  clrs = c("#785DAC","#EC325A","#2AB7CA","#999999")
  
  df.agg2$Condition = factor(df.agg2$Condition,levels=c("HSA.AAG","HSA.Only","AAG.Only","None"))
  
  leg.pos.y = ifelse(grepl("No_Drug",csv.i),.3,.7)
  
  ggplot(df.agg2,aes(x=Protein,y=x.mean,color=Condition,linetype=Condition,size=Condition))+theme_bw()+
    # geom_point(size=2.5)+
    geom_line()+
    geom_errorbar(aes(ymin=x.mean-x.sd,ymax=x.mean+x.sd),size=1.2)+
    xlab("Protein Titration")+
    ylab("Percent Viability")+
    ggtitle(title.i)+
    scale_x_continuous(trans=log2_trans(),
                       breaks=trans_breaks("log2", function(x) 2^x),
                       labels=trans_format("log2", math_format(2^.x)))+
    scale_color_manual(labels=c("HSA+AAG","HSA Only","AAG Only","None"),
                       values=clrs)+
    scale_linetype_manual(values=c("solid","solid","solid","solid"))+
    scale_size_manual(values=c(2.5,2.5,2.5,1.5))+
    theme(
      plot.title = element_text(size=24,face="bold",hjust=0.5),
      axis.title = element_text(size=18,face="bold"),
      axis.text = element_text(size=14,color="black"),
      legend.title = element_text(size=14,face="bold"),
      legend.title.align = 0.5,
      legend.text = element_text(size=12,face="bold"),
      legend.position = c(.2,leg.pos.y),
      legend.background = element_rect(color="black",fill="gray90"),
      legend.key = element_rect(color="black")
    )+
    expand_limits(y=0)+
    guides(linetype=F,size=F)

  # ggheight = ifelse(grepl("No_Drug",csv.i),6,5)
  ggsave(paste(plot.i,".png",sep=""),width=6,height=4,path="Figures")
  
}
