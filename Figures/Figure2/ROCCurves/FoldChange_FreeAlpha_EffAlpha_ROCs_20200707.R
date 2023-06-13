
### ROCs for Binary Classifiers of Resistance ###

## Set up

library(pROC)
library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and clean data

alpha.df = read.csv(file="../../../ReferenceFiles/AlphaValues_True_positive.csv",header=TRUE,stringsAsFactors=FALSE)

# Filter for drug/variant combinations with high confidence (i.e. those with >1 cases)
alpha.df = alpha.df[alpha.df$Actual_Clinical%in%c(0,1),]
alpha.df$Actual_Clinical = as.logical(as.numeric(alpha.df$Actual_Clinical))

# Calculate shift in IC50 over WT
for (drug in unique(alpha.df$drug)) {
  drug.rows = alpha.df$drug == drug
  alpha.df$IC50.shift[drug.rows] = alpha.df$IC50[drug.rows]/alpha.df$IC50[drug.rows&alpha.df$genotype=="WT"]
}

## Build and plot ROC curves
roc.PPB = roc(alpha.df$Actual_Clinical,alpha.df$free.alpha,
              direction=">",levels=c(F,T))
roc.shf = roc(alpha.df$Actual_Clinical,alpha.df$IC50.shift,
              direction="<",levels=c(F,T))
roc.eff = roc(alpha.df$Actual_Clinical,alpha.df$eff.alpha,
              direction=">",levels=c(F,T))

# Plot features
PPB.col = "#ff5a5f"
shf.col = "#f5a300"
eff.col = "#087e8b"

cmmn.theme =   theme(
  plot.title = element_text(size=20,face="bold",hjust=0.5),
  axis.title = element_text(size=18,face="bold"),
  axis.text = element_text(size=14,color="black"),
  legend.title = element_text(size=12,face="bold"),
  legend.title.align = 0.5,
  legend.text = element_text(size=8,face="bold"),
  legend.position = c(.75,.25),
  legend.background = element_rect(color="black",fill="gray95"),
  legend.key = element_rect(color="black")
)

# Plot ROC: IC50 fold-change and free alpha
ggroc(list(roc.shf,roc.PPB),size=2)+theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed",size=2)+
  scale_color_manual("Predictor",labels=c("IC50 Fold Change","Free Alpha"),values=c(shf.col,PPB.col))+
  xlab("Specificity")+ylab("Sensitivity")+
  ggtitle("ROC Curves of\nDrug Kill Rate Classifiers")+
  geom_label(aes(x=0.7,y=0.55,label=paste("AUC =",sprintf("%.2f",round(roc.shf$auc,2)))),size=6,color=shf.col)+
  geom_label(aes(x=0.7,y=0.20,label=paste("AUC =",sprintf("%.2f",round(roc.PPB$auc,2)))),size=6,color=PPB.col)+
  cmmn.theme
  
# ggsave("FreeAlpha_IC50FoldChange_ROCCurve.png",width=5,height=5)

# Plot ROC: IC50 fold-change, free and effective alphas
ggroc(list(roc.eff,roc.shf,roc.PPB),size=2)+theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed",size=2)+
  scale_color_manual("Predictor",labels=c("Effective Alpha","IC50 Fold Change","Free Alpha"),values=c(eff.col,shf.col,PPB.col))+
  xlab("Specificity")+ylab("Sensitivity")+
  ggtitle("ROC Curves of\nDrug Kill Rate Classifiers")+
  geom_label(aes(x=0.85,y=0.85,label=paste("AUC =",sprintf("%.2f",round(roc.eff$auc,2)))),size=6,color=eff.col)+
  geom_label(aes(x=0.67,y=0.57,label=paste("AUC =",sprintf("%.2f",round(roc.shf$auc,2)))),size=6,color=shf.col)+
  geom_label(aes(x=0.25,y=0.5,label=paste("AUC =",sprintf("%.2f",round(roc.PPB$auc,2)))),size=6,color=PPB.col)+
  cmmn.theme

# ggsave("ThreeClassifier_ROCCurve.png",width=5,height=5)
