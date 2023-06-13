
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(reshape2)
library(dr4pl)


## Import and parse data

growth = read.csv('GrowthCurveHeatMap.csv')
subgrowth = growth[,1:12] # Focus analysis on first 13 days, common to all measurements

# Tidy data
df.tidy = melt(subgrowth,id.vars=c("X","X.1"))
colnames(df.tidy) = c("genotype","rep","day","counts")

# Clean data
df.tidy$day = as.numeric(gsub("Day","",df.tidy$day))
df.tidy$genotype = unlist(lapply(strsplit(df.tidy$genotype,"\\("),"[[",1))


## Calculate doubling times

genotypes = unique(df.tidy$genotype)

# Focus only on 1 mL virus replicates for consistency
df.tidy = df.tidy[df.tidy$rep=="1ml virus",]

# Focus on last 3 days (when cells are fully selected)
df.tidy = df.tidy[df.tidy$day>=10,]
df.tidy$day = df.tidy$day-10 # reset day 10 to be day 0

Td.df = data.frame(genotype = genotypes,Td = NA)

for (gen in genotypes) {

  mdl = lm(log(counts)~day,data=df.tidy[df.tidy$genotype==gen,])
  gr = coef(mdl)[2] # growth rate (/day)
  Td = log(2)/gr # doubling time (day)
  
  # Store doubling time
  Td.df$Td[Td.df$genotype==gen] = Td
  
  
}

## Reconstruct IC50 curves and calculate GR50

cmplt.df = read.csv("../../../ReferenceFiles/AlphaValues_True_positive.csv",header=T,stringsAsFactors=F)


# Reorder genotypes
vars = unique(cmplt.df$genotype)
vars = vars[vars!="WT"]
vars = vars[order(substr(vars,2,5))]
cmplt.df$genotype = factor(cmplt.df$genotype,levels=rev(c("WT",vars)))

# Reorder drugs
cmplt.df$drug = factor(cmplt.df$drug,levels=c("Imatinib","Nilotinib","Dasatinib","Bosutinib","Ponatinib"))

doses = 10^seq(-1,5,length.out=100)
T.assay = 3 # duration of assay (days)

cmplt.df$GR50 = NA
cmplt.df$GRhill = NA

for (i in 1:nrow(cmplt.df)) {
  
  genotype = cmplt.df$genotype[i]
  Td.i = Td.df$Td[Td.df$genotype==genotype]
  
  IC50 = cmplt.df$IC50[i]
  hill = cmplt.df$hill[i]
  
  # Reconstruct IC50 Curve
  dose.resp = 1/(1+(doses/IC50)^hill)
  # plot(log10(doses),dose.resp)
  
  # Calculate GR values
  GRs = 2^(1+(log2(dose.resp))/(T.assay/Td.i))-1
  # plot(log10(doses),GRs)

  # Calculate GR50 - note Hafner et al define GR50 where line crosses 0.5, not halfway between GRmin and GRmax
  idx.50 = which(GRs<=0.5)[1]
  idxs = c(idx.50-1,idx.50)
  GR50 = 10^approx(GRs[idxs],log10(doses[idxs]),xout=0.5)$y

  # CalculateHill slope
  fit = dr4pl(doses,GRs)
  GRhill = coef(fit)[3]
  
  cmplt.df$GR50[i] = GR50
  cmplt.df$GRhill[i] = GRhill
  
}

## Calculate GR50 fold change

for (i in 1:nrow(cmplt.df)) {
  
  drug.i = cmplt.df$drug[i]
  
  GR50.WT = cmplt.df$GR50[cmplt.df$drug==drug.i&cmplt.df$genotype=="WT"]
  GR50.var = cmplt.df$GR50[i]
  
  cmplt.df$GR50FC[i] = GR50.var/GR50.WT
  
  IC50.WT = cmplt.df$IC50[cmplt.df$drug==drug.i&cmplt.df$genotype=="WT"]
  IC50.var = cmplt.df$IC50[i]
  
  cmplt.df$IC50FC[i] = IC50.var/IC50.WT
  
}

## Plot heat map

# Common theme
cmmn.thm = theme(
  plot.title = element_text(size=20,face="bold",hjust=0.5),
  axis.title = element_text(size=18,face="bold"),
  axis.text = element_text(size=12,color="black",face="bold"),
  axis.text.x = element_text(angle=15,vjust=0.5),
  legend.title = element_text(size=14,face="bold"),
  legend.title.align = 0.5,
  legend.text = element_text(size=14,face="bold")
)

# GR50 fold change
ggplot(cmplt.df,aes(x=drug,y=genotype,fill=GR50FC))+theme_bw()+
  geom_raster()+
  xlab("Drugs")+ylab("Genotype")+
  ggtitle("Fold Change in GR50 over WT")+
  scale_fill_gradientn("GR50\nFold-Change",colors=rev(brewer.pal(10,"RdYlBu")),
                       trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x)))+
  cmmn.thm
# ggsave("GR50Heatmap.png",width=6,height=5)

## Build classifiers and ROC curves

# Filter for drug/variant combinations with high confidence (i.e. those with >1 cases)
cmplt.df = cmplt.df[cmplt.df$Actual_Clinical%in%c(0,1),]
cmplt.df$Actual_Clinical = as.logical(as.numeric(cmplt.df$Actual_Clinical))


roc.IC50 = roc(cmplt.df$Actual_Clinical,cmplt.df$IC50FC,
               direction="<",levels=c(F,T))
roc.GR50 = roc(cmplt.df$Actual_Clinical,cmplt.df$GR50FC,
           direction="<",levels=c(F,T))

# Plot features
IC50.col = "#f5a300"
GR50.col = "purple"

ggroc(list(roc.IC50,roc.GR50),size=2)+theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed",size=2)+
  scale_color_manual("Predictor",labels=c("IC50 Fold Change","GR50 Fold Change"),values=c(IC50.col,GR50.col))+
  xlab("Specificity")+ylab("Sensitivity")+
  ggtitle("ROC Curves of\nDrug Kill Rate Classifiers")+
  geom_label(aes(x=0.7,y=0.55,label=paste("AUC =",sprintf("%.2f",round(roc.IC50$auc,2)))),size=6,color=IC50.col)+
  geom_label(aes(x=0.7,y=0.20,label=paste("AUC =",sprintf("%.2f",round(roc.GR50$auc,2)))),size=6,color=GR50.col)+
  cmmn.thm
# ggsave("GR50_IC50_ROC.png",width=7,height=5)

