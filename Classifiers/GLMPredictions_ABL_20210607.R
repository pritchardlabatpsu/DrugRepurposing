
## Predicting Resistance from Logistic Regression Models ##

## Set up

library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load and clean files

# Load in models (see ResistanceLogisticRegression file)
mdls = readRDS("LogisticRegressionModels.rds")

# Load in data
alpha.df = read.csv("AlphaValues_True_positive.csv",header=T,stringsAsFactors=F)

# Filter for drug/variant combinations with high confidence (i.e. those with >1 cases)
alpha.df$Actual_Clinical[alpha.df$Actual_Clinical==0] = F
alpha.df$Actual_Clinical[alpha.df$Actual_Clinical==1] = T
alpha.df$Actual_Clinical[alpha.df$Actual_Clinical=="N"] = NA

## Predict resistance

# Classification from best model (eff.alpha chosen based on AIC, BIC, LOOCV and for parsimony)
mdl = mdls$eff.alpha
alpha.df$pred.response = predict(mdl,
                                 newdata=list(eff.alpha=alpha.df$eff.alpha),
                                 type="response")
alpha.df$pred.clinical = alpha.df$pred.response>=0.5
alpha.df$correct.pred = alpha.df$Actual_Clinical==alpha.df$pred.clinical

# Find threshold for eff.alpha that gives probability of resistance = 0.5 (cutoff value)
eff.alpha.vec = seq(min(alpha.df$eff.alpha),max(alpha.df$eff.alpha),by=0.01)
log.pred = predict(mdl,
                   newdata=list(eff.alpha=eff.alpha.vec),
                   type="response")
thrshld = eff.alpha.vec[log.pred<=0.5][1] # threshold at eff.alpha = 0.44

# What is the associated imatinib concentration for this kill rate?
t.assay = 3 # duration of assay [days]
y.thrshld = exp(-thrshld*t.assay)

IC50 = alpha.df$IC50[alpha.df$drug=="Imatinib"&alpha.df$genotype=="WT"]
hill = alpha.df$hill[alpha.df$drug=="Imatinib"&alpha.df$genotype=="WT"]

dose.thrshld = IC50*(1/y.thrshld-1)^(1/hill) # effective dose must be >= 281 nM to target WT ABL

# What are the TPR and TNR?

sub.df = alpha.df[!is.na(alpha.df$Actual_Clinical),] # subset on drug/variant combos where clinical data are available
sub.df$Actual_Clinical = as.logical(sub.df$Actual_Clinical)
TPR = sum(sub.df$Actual_Clinical&sub.df$pred.clinical)/sum(sub.df$Actual_Clinical) # 90.3% sensitive
TNR = sum(!sub.df$Actual_Clinical&!sub.df$pred.clinical)/sum(!sub.df$Actual_Clinical) # 84.2% specific

## Plot results

# Order drugs and genotypes
alpha.df$drug = factor(alpha.df$drug,levels=c("Imatinib","Nilotinib","Dasatinib","Bosutinib","Ponatinib"))

vars = unique(alpha.df$genotype)
vars = vars[vars!="WT"]
vars = vars[order(substr(vars,2,5))]
alpha.df$genotype = factor(alpha.df$genotype,levels=rev(c("WT",vars)))

alpha.df$correct.pred[is.na(alpha.df$correct.pred)] = "NA"
alpha.df$correct.pred = factor(alpha.df$correct.pred,levels=c("TRUE","FALSE","NA"))

alpha.df$resistant.str = ifelse(alpha.df$Actual_Clinical,"(R)","(S)")

ggplot(alpha.df,aes(x=drug,y=genotype,fill=correct.pred))+theme_minimal()+
  geom_tile(color="black",size=1.25)+
  geom_text(aes(label=paste(resistant.str,sprintf("%.2f",eff.alpha))),size=5)+
  xlab("Drug")+ylab("Genotype")+
  ggtitle("Predicted Clinical Outcomes: ABL")+
  scale_fill_manual("Predicted",
                    labels = c("Correct","Incorrect","NA"),
                    values=c("#00a7e1","#ff4b33","#dedede"))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=14,color="black",face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size=14,face="bold"),
    legend.position = "bottom"
  )
  
# ggsave("ABLHeatMap.png",height=12,width=8)

### Predict Resistance using free Cave

# Clear environment
rm(list=ls())

# Load and clean files

# Load in models (see ResistanceLogisticRegression file)
mdls = readRDS("LogisticRegressionModels.rds")

# Load in data
alpha.df = read.csv("AlphaValues_True_positive.csv",header=T,stringsAsFactors=F)

# Filter for drug/variant combinations with high confidence (i.e. those with >1 cases)
alpha.df$Actual_Clinical[alpha.df$Actual_Clinical==0] = F
alpha.df$Actual_Clinical[alpha.df$Actual_Clinical==1] = T
alpha.df$Actual_Clinical[alpha.df$Actual_Clinical=="N"] = NA

# Classification from free Cave model
mdl = mdls$free.alpha
alpha.df$pred.response = predict(mdl,
                                 newdata=list(free.alpha=alpha.df$free.alpha),
                                 type="response")
alpha.df$pred.clinical = alpha.df$pred.response>=0.5
alpha.df$correct.pred = alpha.df$Actual_Clinical==alpha.df$pred.clinical

# Find threshold for eff.alpha that gives probability of resistance = 0.5 (cutoff value)
free.alpha.vec = seq(min(alpha.df$eff.alpha),max(alpha.df$eff.alpha),by=0.01)
log.pred = predict(mdl,
                   newdata=list(free.alpha=free.alpha.vec),
                   type="response")
thrshld = free.alpha.vec[log.pred<=0.5][1] # threshold at free.alpha near 0

# What is the associated imatinib concentration for this kill rate?
t.assay = 3 # duration of assay [days]
y.thrshld = exp(-thrshld*t.assay)

IC50 = alpha.df$IC50[alpha.df$drug=="Imatinib"&alpha.df$genotype=="WT"]
hill = alpha.df$hill[alpha.df$drug=="Imatinib"&alpha.df$genotype=="WT"]

dose.thrshld = IC50*(1/y.thrshld-1)^(1/hill) # effective dose must be >= 281 nM to target WT ABL

# What are the TPR and TNR?

sub.df = alpha.df[!is.na(alpha.df$Actual_Clinical),] # subset on drug/variant combos where clinical data are available
sub.df$Actual_Clinical = as.logical(sub.df$Actual_Clinical)
TPR = sum(sub.df$Actual_Clinical&sub.df$pred.clinical)/sum(sub.df$Actual_Clinical) # 90.3% sensitive
TNR = sum(!sub.df$Actual_Clinical&!sub.df$pred.clinical)/sum(!sub.df$Actual_Clinical) # 84.2% specific

## Plot results

# Order drugs and genotypes
alpha.df$drug = factor(alpha.df$drug,levels=c("Imatinib","Nilotinib","Dasatinib","Bosutinib","Ponatinib"))

vars = unique(alpha.df$genotype)
vars = vars[vars!="WT"]
vars = vars[order(substr(vars,2,5))]
alpha.df$genotype = factor(alpha.df$genotype,levels=rev(c("WT",vars)))

alpha.df$correct.pred[is.na(alpha.df$correct.pred)] = "NA"
alpha.df$correct.pred = factor(alpha.df$correct.pred,levels=c("TRUE","FALSE","NA"))

alpha.df$resistant.str = ifelse(alpha.df$Actual_Clinical,"(R)","(S)")

ggplot(alpha.df,aes(x=drug,y=genotype,fill=correct.pred))+theme_minimal()+
  geom_tile(color="black",size=1.25)+
  geom_text(aes(label=paste(resistant.str,sprintf("%.2f",eff.alpha))),size=5)+
  xlab("Drug")+ylab("Genotype")+
  ggtitle("Free Cave Predicted Clinical Outcomes")+
  scale_fill_manual("Predicted",
                    labels = c("Correct","Incorrect","NA"),
                    values=c("#00a7e1","#ff4b33","#dedede"))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=14,color="black",face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size=14,face="bold"),
    legend.position = "bottom"
  )

# ggsave("ABLHeatMap_Free.png",height=12,width=8)
