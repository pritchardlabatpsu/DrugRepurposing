
### Predicting resistance in EGFR and KIT ###

## Set up

library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Loop through files and plot heat maps

path = "./"

files = c("AlphaValues_EGFRCompiled_20210729.csv","AlphaValues_KITCompiled_081920.csv")

# Logistic models (see Figure 3)
mdls = readRDS("LogisticRegressionModels.rds")
predictors = c("eff.alpha","free.alpha")

## Read in and plot data

# Dimensions for saving plots
dim = list(width=c(7,6),height=c(5.5,7))

cmplt.df = data.frame()

for (i in 1:length(files)) {
  
  # Read in data
  file.i = paste(path,files[i],sep="")
  alpha.df = read.csv(file.i,header=T,stringsAsFactors=F)
  GOI = strsplit(file.i,"_")[[1]][2]
  GOI = gsub("Compiled","",GOI)
  
  for (predictor in predictors) {
    
    mdl = mdls[predictor][[1]]
  
    # Predict resistance using logistic model
    alpha.df$pred.response = predict(mdl,
                                     newdata=list(eff.alpha=alpha.df$eff.Cave.alpha,
                                                  free.alpha=alpha.df$free.Cave.alpha),
                                     type="response")
    alpha.df$pred.clinical = alpha.df$pred.response>=0.5
    alpha.df$correct.pred = alpha.df$Actual_Clinical==alpha.df$pred.clinical
    
    alpha.df$genotype = paste(alpha.df$primary,alpha.df$secondary)
    alpha.df$genotype = factor(alpha.df$genotype,levels = rev(unique(alpha.df$genotype)))
    
    alpha.df$resistant.str = ifelse(alpha.df$Actual_Clinical,"(R)","(S)")
    alpha.df$resistant.str[is.na(alpha.df$resistant.str)] = "(NA)"
    
    if (predictor=="eff.alpha") {
      cmplt.df = rbind(cmplt.df,alpha.df)
    }
    
    ## Plot results
    
    # Fill in NAs
    expnd.df = expand.grid(drug=unique(alpha.df$drug),genotype=unique(alpha.df$genotype))
    for (j in 1:nrow(alpha.df)) {
      drug.j = alpha.df$drug[j]
      genotype.j = alpha.df$genotype[j]
      expnd.df$eff.Cave.alpha[expnd.df$drug==drug.j&expnd.df$genotype==genotype.j] = alpha.df$eff.Cave.alpha[j]
      expnd.df$free.Cave.alpha[expnd.df$drug==drug.j&expnd.df$genotype==genotype.j] = alpha.df$free.Cave.alpha[j]
      expnd.df$correct.pred[expnd.df$drug==drug.j&expnd.df$genotype==genotype.j] = alpha.df$correct.pred[j]
      expnd.df$resistant.str[expnd.df$drug==drug.j&expnd.df$genotype==genotype.j] = alpha.df$resistant.str[j]
      
    }
    
    expnd.df$correct.pred[is.na(expnd.df$correct.pred)] = "NA"
    expnd.df$correct.pred = factor(expnd.df$correct.pred,levels=c("TRUE","FALSE","NA"))
    
    # Plot alpha values
    if (predictor=="eff.alpha") {
      plt.alpha = expnd.df$eff.Cave.alpha
    } else {
      plt.alpha = expnd.df$free.Cave.alpha
    }
    
    sfm = scale_fill_manual("Predicted",
                            breaks = c("TRUE","FALSE","NA"),
                            labels = c("Correct","Incorrect","NA"),
                            values=c("#00a7e1","#ff4b33","#dedede"))
    
    ggplot(expnd.df,aes(x=drug,y=genotype,fill=correct.pred))+theme_minimal()+
      geom_tile(color="black",size=1.25)+
      geom_text(aes(label=paste(resistant.str,sprintf("%.2f",plt.alpha))),size=5)+
      xlab("Drug")+ylab("Genotype\n")+
      ggtitle(paste("Predicted Clinical Outcomes:\n",GOI))+
      sfm+
      theme(
        plot.title = element_text(size=20,face="bold",hjust=0.5),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=14,color="black",face="bold"),
        axis.text.y = element_text(hjust=0),
        legend.title = element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.position = "bottom"
      )
    
    # ggsave(paste(GOI,gsub("\\.","",predictor),"HeatMap.png",sep=""),width=dim$width[i],height=dim$height[i])

  }
}

## Sensitivity/Specificity of Combined Sets

sens = sum(cmplt.df$Actual_Clinical&cmplt.df$correct.pred,na.rm=T)/sum(cmplt.df$Actual_Clinical,na.rm=T) # (number of correctly classified resistance cases)/(total number of resistance cases)
spec = sum(!cmplt.df$Actual_Clinical&cmplt.df$correct.pred,na.rm=T)/sum(!cmplt.df$Actual_Clinical,na.rm=T) # (number of correctly classified sensitive cases)/(total number of sensitive cases)
