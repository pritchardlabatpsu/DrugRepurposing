

### ROCs for Binary Classifiers of Resistance ###

## Set up

library(pROC)
library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and clean data

alpha.df = read.csv(file= "../../../ReferenceFiles/AlphaValues_True_positive.csv",header=TRUE,stringsAsFactors=FALSE)

# Reorder drugs
alpha.df$drug = factor(alpha.df$drug,levels=c("Imatinib","Nilotinib","Dasatinib","Bosutinib","Ponatinib"))

# Reorder genotypes
vars = unique(alpha.df$genotype)
vars = vars[vars!="WT"]
vars = vars[order(substr(vars,2,5))]
alpha.df$genotype = factor(alpha.df$genotype,levels=rev(c("WT",vars)))

## Calculate shift in IC50 over WT
for (drug in unique(alpha.df$drug)) {
  drug.rows = alpha.df$drug == drug
  alpha.df$IC50.shift[drug.rows] = alpha.df$IC50[drug.rows]/alpha.df$IC50[drug.rows&alpha.df$genotype=="WT"]
}

## Calculate drug kill rate (alpha) at Cave (uncorrected)

# Note Cave in nM (uncorrected - Figure 6c)
alpha.df$Cave = NA
alpha.df$Cave[alpha.df$drug=="Imatinib"] = 3377
alpha.df$Cave[alpha.df$drug=="Nilotinib"] = 2754
alpha.df$Cave[alpha.df$drug=="Dasatinib"] = 27
alpha.df$Cave[alpha.df$drug=="Bosutinib"] = 287
alpha.df$Cave[alpha.df$drug=="Ponatinib"] = 101

alpha.df$Cave.alpha = 1/(1+(alpha.df$Cave/alpha.df$IC50)^alpha.df$hill)

## Plot heatmaps

alpha.df$Actual_Clinical = as.logical(as.numeric(alpha.df$Actual_Clinical))

## Build and plot ROC curves
roc.PPB = roc(alpha.df$Actual_Clinical,alpha.df$free.alpha,
              direction=">",levels=c(F,T))
roc.shf = roc(alpha.df$Actual_Clinical,alpha.df$IC50.shift,
              direction="<",levels=c(F,T))
roc.Cave = roc(alpha.df$Actual_Clinical,alpha.df$Cave.alpha,
               direction=">",levels=c(F,T))

# Plot features
PPB.col = "#ff5a5f"
Cave.col = "#13ab02"
  
cmmn.thm = theme(
  plot.title = element_text(size=20,face="bold",hjust=0.5),
  axis.title = element_text(size=18,face="bold"),
  axis.text = element_text(size=12,color="black",face="bold"),
  axis.text.x = element_text(angle=15,vjust=0.5),
  legend.title = element_text(size=14,face="bold"),
  legend.title.align = 0.5,
  legend.text = element_text(size=14,face="bold")
)

## Plot heatmap

ggplot(alpha.df,aes(x=drug,y=genotype,fill=Cave.alpha))+theme_bw()+
  geom_raster()+
  xlab("Drugs")+ylab("Genotype")+
  ggtitle("Drug Kill Rate: Uncorrected Alpha")+
  scale_fill_gradientn("Uncorrected\nAlpha",colors=rev(brewer.pal(10,"RdYlBu")))+
  cmmn.thm

# ggsave("UncorrectedAlpha_Heatmap.png",width=6,height=5)

## Plot ROC: free Cave and uncorrected Cave
ggroc(list(roc.PPB,roc.Cave),size=2)+theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed",size=2)+
  scale_color_manual("Predictor",labels=c("Free Alpha","Uncorrected Alpha"),values=c(PPB.col,Cave.col))+
  xlab("Specificity")+ylab("Sensitivity")+
  ggtitle("ROC Curves of\nDrug Kill Rate Classifiers")+
  geom_label(aes(x=0.67,y=0.57,label=paste("AUC =",sprintf("%.2f",round(roc.PPB$auc,2)))),size=6,color=PPB.col)+
  geom_label(aes(x=0.25,y=0.5,label=paste("AUC =",sprintf("%.2f",round(roc.Cave$auc,2)))),size=6,color=Cave.col)+
  cmmn.thm

# ggsave("FreeAlpha_UncorrectedAlpha_ROC.png",width=7,height=5)
