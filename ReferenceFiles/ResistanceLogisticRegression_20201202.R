
### Logistic Regression Models of Resistance Classification ###

## Set up

library(stringr)
library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and clean data

alpha.df = read.csv("AlphaValues_True_positive.csv",header=T,stringsAsFactors=F)

# Filter for drug/variant combinations with high confidence (i.e. those with >1 cases)
alpha.df = alpha.df[alpha.df$Actual_Clinical%in%c(0,1),]
alpha.df$Actual_Clinical = as.logical(as.numeric(alpha.df$Actual_Clinical))

# Volumes of distribution
Vd = c(252,1223,2500,800,6080)
names(Vd) = c("Imatinib","Ponatinib","Dasatinib","Nilotinib","Bosutinib")

# Fold change in IC50 over WT
for (drug in unique(alpha.df$drug)) {
  
  alpha.df$IC50.fold.change[alpha.df$drug==drug] = alpha.df$IC50[alpha.df$drug==drug]/alpha.df$IC50[alpha.df$drug==drug&alpha.df$genotype=="WT"]
  
}

for (drug in names(Vd)) {
  alpha.df$Vd[alpha.df$drug==drug] = Vd[drug]
}

## Build logistic regression models

# Define variables used in each model
all.vars = c("free.alpha","eff.alpha","eff.max.alpha","Vd")

# Develop all combinations of variables
vars.comb = c()
r = 1
for (i in 1:length(all.vars)) {
  comb.mat = combn(all.vars,i)
  for (j in 1:ncol(comb.mat)) {
    vars.comb[r] = paste(comb.mat[,j],collapse="+")
    r = r+1
  }
}

# Build and evaluate logistic regression models
mdl.smmry = as.data.frame(matrix(nrow=length(vars.comb),ncol=5))
mdls = list()
colnames(mdl.smmry) = c("vars","nvars","AIC","BIC","LOOCV.err")

for (i in 1:length(vars.comb)) {
  
  vars.i = vars.comb[i]
  mdl.smmry$vars[i] = vars.i
  mdl.smmry$nvars[i] = str_count(vars.i,"\\+")+1
  
  frmla = as.formula(paste("Actual_Clinical ~",vars.i))
  mdl = glm(formula = frmla,
            data = alpha.df,
            family = binomial)
  mdls[[vars.i]] = mdl
  
  mdl.smmry$AIC[i] = AIC(mdl)
  mdl.smmry$BIC[i] = BIC(mdl)
  
  # Calculate LOOCV error
  CV.pred = rep(NA,nrow(alpha.df))
  
  for (j in 1:nrow(alpha.df)) {
    mdl.CV = glm(formula = frmla,
                 data = alpha.df[-j,],
                 family = binomial)
    pred.j = predict(mdl.CV,
                     newdata = alpha.df[j,],
                     type="response")
    CV.pred[j] = pred.j>=0.5
  }
  
  mdl.smmry$LOOCV.err[i] = mean(CV.pred!=alpha.df$Actual_Clinical)

}

# saveRDS(mdls,file="LogisticRegressionModels.rds")

# Find minimum criterion for each metric
mdl.smmry$AIC.min = F
mdl.smmry$BIC.min = F
mdl.smmry$LOOCV.err.min = F

for (i in 1:length(all.vars)) {
  mdl.smmry$AIC.min[mdl.smmry$nvars==i&
                      mdl.smmry$AIC==min(mdl.smmry$AIC[mdl.smmry$nvars==i])] = T
  mdl.smmry$BIC.min[mdl.smmry$nvars==i&
                      mdl.smmry$BIC==min(mdl.smmry$BIC[mdl.smmry$nvars==i])] = T
  mdl.smmry$LOOCV.err.min[mdl.smmry$nvars==i&
                          mdl.smmry$LOOCV.err==min(mdl.smmry$LOOCV.err[mdl.smmry$nvars==i])] = T
}

# # Visualize model fits
# plot(alpha.df$eff.alpha,alpha.df$Actual_Clinical)
# yfit = predict(mdls$`free.alpha+Vd`,type="response")
# ordr = order(alpha.df$eff.alpha)
# points(alpha.df$eff.alpha[ordr],yfit[ordr])

## Plot results

cmmn.theme =   theme(
  plot.title = element_text(size=20,face="bold",hjust=0.5),
  axis.title = element_text(size=18,face="bold"),
  axis.text = element_text(size=14,color="black")
)

# Plot AIC
AIC.col = "#3c5190"
# ggplot(mdl.smmry[mdl.smmry$AIC.min,],aes(x=nvars,y=AIC))+theme_bw()+
#   geom_point(size=5,color=AIC.col)+
#   geom_line(size=3,color=AIC.col)+
#   xlab("Best N-Variable Model")+
#   ggtitle("Model Selection: AIC")+
#   cmmn.theme

# alt:
ggplot(mdl.smmry,aes(x=nvars,y=AIC))+theme_bw()+
  geom_point(aes(color=AIC.min,alpha=AIC.min),size=5)+
  geom_line(data=mdl.smmry[mdl.smmry$AIC.min,],size=3,color=AIC.col)+
  xlab("N-Variable Models")+
  ggtitle("Model Selection: AIC")+
  scale_color_manual(values=c("gray",AIC.col))+
  scale_alpha_manual(values=c(0.75,1))+
  guides(color=F,alpha=F)+
  cmmn.theme

# ggsave("ModelSelection_AIC.png",width=6,height=4)

# Plot BIC
BIC.col = "#496a81"
# ggplot(mdl.smmry[mdl.smmry$BIC.min,],aes(x=nvars,y=BIC))+theme_bw()+
#   geom_point(size=5,color=BIC.col)+
#   geom_line(size=3,color=BIC.col)+
#   xlab("Best N-Variable Model")+
#   ggtitle("Model Selection: BIC")+
#   cmmn.theme

# alt:
ggplot(mdl.smmry,aes(x=nvars,y=BIC))+theme_bw()+
  geom_point(aes(color=BIC.min,alpha=BIC.min),size=5)+
  geom_line(data=mdl.smmry[mdl.smmry$BIC.min,],size=3,color=BIC.col)+
  xlab("N-Variable Models")+
  ggtitle("Model Selection: BIC")+
  scale_color_manual(values=c("gray",BIC.col))+
  scale_alpha_manual(values=c(0.75,1))+
  guides(color=F,alpha=F)+
  cmmn.theme

# ggsave("ModelSelection_BIC.png",width=6,height=4)

# Plot LOOCV
CV.col = "#66999b"
# ggplot(mdl.smmry[mdl.smmry$LOOCV.err.min,],aes(x=nvars,y=LOOCV.err))+theme_bw()+
#   geom_point(size=5,color=CV.col)+
#   geom_line(size=3,color=CV.col)+
#   xlab("Best N-Variable Model")+
#   ylab("Cross Validation Error")+
#   ggtitle("Model Selection: LOOCV")+
#   cmmn.theme

# alt:
ggplot(mdl.smmry,aes(x=nvars,y=LOOCV.err))+theme_bw()+
  geom_point(aes(color=LOOCV.err.min,alpha=LOOCV.err.min),size=5)+
  geom_line(data=mdl.smmry[mdl.smmry$LOOCV.err.min,],size=3,color=CV.col)+
  xlab("N-Variable Models")+
  ylab("Cross Validation Error")+
  ggtitle("Model Selection: LOOCV")+
  scale_color_manual(values=c("gray",CV.col))+
  scale_alpha_manual(values=c(0.75,1))+
  guides(color=F,alpha=F)+
  cmmn.theme

# ggsave("ModelSelection_LOOCV.png",width=6,height=4)

