
### Plot prospective predicitons for avaprtinib and ripretinib ###

## Set up

library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

## Import data

KIT.df = read.csv("ProspectiveKITCompiled_090320.csv")

## Aggregate data

# KIT.df$secondary[KIT.df$secondary=="WT"] = ""

KIT.df$genotype = paste(KIT.df$primary,KIT.df$secondary)
KIT.df$genotype[KIT.df$primary=="V560D"&KIT.df$secondary=="V560D"] = "V560D WT"

## Make predictions (model from Fig 3)

mdls = readRDS("LogisticRegressionModels.rds")
mdl = mdls$eff.alpha

KIT.df$pred.response = predict(mdl,
                               newdata=list(eff.alpha=KIT.df$eff.Cave.alpha),
                               type="response")
KIT.df$pred.clinical = KIT.df$pred.response>=0.5

## Plot results

# Order aesthetics
KIT.df$pred.clinical = factor(KIT.df$pred.clinical,levels=c(T,F))

gen.order = c("Ex9Mut WT","Ex9Mut V654A","Ex9Mut D816H",
              "Ex11Mut WT","Ex11Mut V560D","Ex11Mut V654A","Ex11Mut T670I","Ex11Mut D816G","Ex11Mut D816H","Ex11Mut D816V","Ex11Mut D816Y","Ex11Mut D820A","Ex11Mut D820G","Ex11Mut N822K","Ex11Mut Y823D","Ex11Mut A829P",
              "V560D WT","V560D V654A","V560D D816H")
KIT.df$genotype = factor(KIT.df$genotype,
                         levels = rev(gen.order))

# Plotting
ggplot(KIT.df,aes(x=drug,y=genotype,fill=pred.clinical))+theme_minimal()+
  geom_tile(color="black",size=1.25)+
  geom_text(aes(label=sprintf("%.2f",eff.Cave.alpha),color=pred.clinical),size=5)+
  xlab("Drug")+ylab("Genotype\n")+
  ggtitle(paste("Predicted Clinical Outcomes:\nKIT"))+
  scale_fill_manual("Predicted Phenotype: ",
                    breaks = c("FALSE","TRUE"),
                    labels = c("Sensitive","Resistant"),
                    values=c("#49216E","#BFA8C7"))+
  scale_color_manual(values=c("white","black"))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=14,color="black",face="bold"),
    axis.text.y = element_text(hjust=0),
    legend.title = element_text(size=18,face="bold",hjust=0.5),
    legend.text = element_text(size=14,face="bold"),
    legend.position = "bottom"
  )+
  guides(color=F)

ggsave("ProspectiveKITPredictionHeatMap.png",width=6,height=5)

## Subset on Avapritinib/Ex11

ava11.df = KIT.df[KIT.df$drug=="Avapritinib"&KIT.df$primary=="Ex11Mut",]

# Reorder mutations
scndry = unique(ava11.df$secondary)
vars = as.character(scndry[scndry!="WT"])
ordr = order(as.numeric(substr(vars,2,4)))
scndry.ordrd = c("WT",vars[ordr])

ava11.df$secondary = factor(ava11.df$secondary,levels=rev(scndry.ordrd))

ggplot(ava11.df,aes(x=drug,y=secondary,fill=pred.clinical))+theme_minimal()+
  geom_tile(color="black",size=1.25)+
  geom_text(aes(label=sprintf("%.2f",eff.Cave.alpha),color=pred.clinical),size=5)+
  xlab("Drug")+ylab("Genotype\n")+
  ggtitle(paste("Predicted Clinical Outcomes:\nKIT Ex11Mut"))+
  scale_fill_manual("Predicted\nPhenotype: ",
                    breaks = c("FALSE","TRUE"),
                    labels = c("Sensitive","Resistant"),
                    values=c("#4E707E","#C0D1D8"))+
  scale_color_manual(values=c("white","black"))+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=14,color="black",face="bold"),
    axis.text.y = element_text(hjust=0),
    legend.title = element_text(size=18,face="bold",hjust=0.5),
    legend.title.align = 0.5,
    legend.text = element_text(size=14,face="bold",hjust=0.5)
  )+
  guides(color=F)

# ggsave("ProspectiveKIT11PredictionHeatMap.png",width=5,height=5)
