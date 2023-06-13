
### Plot Imatinib Preclinical/Clinical Summary ###

## Set up

library(ggplot2)
library(scales)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import data

df = read.csv("../../ReferenceFiles/RepurposingClinicalTrials/ImatinibClinicalTrials_20210214.csv",header=T)

## Drop repeated diseases (consider highest dose evaluated)

df = df[df$Simplified!="Repeat Disease",]

df$Disease = factor(df$Disease,levels=df$Disease[order(df$Lowest.Dose.uM)])
df$Simplified = factor(df$Simplified,levels=df$Simplified[order(df$Lowest.Dose.uM)])

df$Approved[df$Approved=="TRUE"] = "Approved"
df$Approved[df$Approved=="FALSE"] = "Not Approved"
df$Approved = factor(df$Approved,levels=c("Approved","Not Approved"))

## Plot results

ggplot(df,aes(x=Disease,y=Lowest.Dose.uM,color=Status))+theme_bw()+
  geom_point(size=10)+
  facet_grid(~Approved,scales="free_x",space="free_x")+
  xlab("Unique Diseases")+ylab("Lowest Efficacious Dose\nin in vitro System [uM]")+
  ggtitle("Imatinib Repurposing:\nPreclinical Efficacy and Clinical Outcomes")+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)))+
  scale_color_manual("Clinical\nStatus",values=c("#F53D00","#FDC835","#33673B"))+
  theme(
    plot.title = element_text(size=38,color="black",face="bold",hjust=0.5),
    axis.title = element_text(size=30,color="black",face="bold"),
    axis.text = element_text(color="black",face="bold"),
    axis.text.y = element_text(size=24),
    axis.text.x = element_text(size=16,angle=75,vjust=1,hjust=1),
    legend.title = element_text(size=30,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=24,face="bold"),
    strip.text = element_text(size=30,face="bold"),
    strip.background = element_rect(fill="gray95")
  )

# ggsave("ImatinibPreclinicalClinicalSummary.png",width=18,height=18)

# Alternatively, plot so scale matches range of histogram of preclinical 2018 data

ggplot(df,aes(x=Disease,y=Lowest.Dose.uM,size=Enrollment,color=Status))+theme_bw()+
  geom_hline(yintercept = .444,size=1.5,color="#087E8B",linetype="dashed")+
  geom_hline(yintercept = .579,size=1.5,color="#044c54",linetype="dashed")+
  geom_point()+
  # geom_label(aes(label=Phase,fill=Status))+
  facet_grid(~Approved,scales="free_x",space="free_x")+
  xlab("Unique Diseases")+ylab("Lowest Efficacious Dose\nin in vitro System [uM]")+
  ggtitle("Imatinib Repurposing:\nPreclinical Efficacy and Clinical Outcomes")+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)),limits = c(10^(-2.5),10^(2.5)))+
  scale_color_manual("Clinical\nStatus",values=c("#F53D00","#FDC835","#33673B"))+
  # scale_color_manual("Clinical\nStatus",values=c("white","black","white"))+
  scale_size_continuous(range = c(6,15))+
  theme(
    plot.title = element_text(size=38,color="black",face="bold",hjust=0.5),
    axis.title = element_text(size=30,color="black",face="bold"),
    axis.text = element_text(color="black",face="bold"),
    axis.text.y = element_text(size=24),
    axis.text.x = element_text(size=16,angle=75,vjust=1,hjust=1),
    # axis.text.x = element_blank(),
    legend.title = element_text(size=30,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=24,face="bold"),
    strip.text = element_text(size=30,face="bold"),
    strip.background = element_rect(fill="gray95")
  )+
  guides(color=guide_legend(override.aes=list(size=12)))

# ggsave("ImatinibPreclinicalClinicalSummaryAlt.png",width=18,height=18)

## Alternatively, segregate by max dose used in evaluation

df.dose = df[df$Max.Dose!="340 mg/m2 QD",] # Remove pediatric PVS, since dosing is based on surface area

df.dose$Max.Dose = factor(df.dose$Max.Dose,levels=c("300 mg QD","400 mg QD","600 mg QD","400 mg BID"))

imat.dose = read.csv("../../ReferenceFiles/PreclinicalMetaAnalysis/ImatinibDosing_20221011.csv",header=T)

for (i in 1:nrow(df.dose)) {
  df.dose$eff.Cave.uM[i] = imat.dose$eff.Cave.nM[imat.dose$Dose==df.dose$Max.Dose[i]]/1000
  df.dose$eff.Cmax.uM[i] = imat.dose$eff.Cmax.nM[imat.dose$Dose==df.dose$Max.Dose[i]]/1000
}

g = ggplot(df.dose,aes(x=Simplified,y=Lowest.Dose.uM,size=Enrollment,color=Status))+theme_bw()+
  facet_grid(~Max.Dose,scales="free_x",space="free_x")+
  geom_hline(aes(yintercept=eff.Cave.uM),size=1.5,color="#087E8B",linetype="dashed")+
  geom_hline(aes(yintercept=eff.Cmax.uM),size=1.5,color="#04464D",linetype="dashed")+
  geom_point()+
  # geom_label(aes(label=Phase,fill=Status))+
  xlab("Unique Diseases")+ylab("Lowest Efficacious Dose\nin in vitro System [uM]")+
  ggtitle("Imatinib Repurposing:\nPreclinical Efficacy and Clinical Outcomes")+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x)),limits = c(10^(-2.5),10^(2.5)))+
  scale_color_manual("Clinical\nStatus",values=c("#F53D00","#FDC835","#33673B"))+
  # scale_color_manual("Clinical\nStatus",values=c("white","black","white"))+
  scale_size_continuous(range = c(6,15))+
  theme(
    plot.title = element_text(size=38,color="black",face="bold",hjust=0.5),
    axis.title = element_text(size=30,color="black",face="bold"),
    axis.text = element_text(color="black",face="bold"),
    axis.text.y = element_text(size=24),
    axis.text.x = element_text(size=16,angle=45,vjust=1,hjust=1),
    # axis.text.x = element_blank(),
    legend.title = element_text(size=30,face="bold"),
    legend.title.align = 0.5,
    legend.text = element_text(size=24,face="bold"),
    strip.text = element_text(size=30,face="bold"),
    strip.background = element_rect(fill="gray95")
  )+
  guides(color=guide_legend(override.aes=list(size=12)))

# Adjust facet width
gp = ggplotGrob(g)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
gp$widths[facet.columns][[1]] = gp$widths[facet.columns][[1]] * 4 # increase width of 300 mg QD panel by factor of 3

# Plot result
grid::grid.draw(gp)

# pdf("ImatinibPreclinicalClinicalSummaryByDose_20221011.pdf",height=18,width=20)
# grid::grid.draw(gp)
# dev.off()
