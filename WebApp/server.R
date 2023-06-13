library(shiny)
library(ggplot2)
library(dr4pl)

### Define server logic
server = function(input, output) {
  
  ## Reactive expression
  dr_df <- eventReactive( input$submit_button, {
    read.table(input$dr_data$datapath,
               header = F, sep = ",", stringsAsFactors = F) 
  })
  
  ## Plot output
  output$dr_plot = renderPlot({
    
    # Read in data
    df = dr_df()
    colnames(df) = c("conc","serum","replicate","rel.viab")
    
    # Transform drug concentrations
    if(input$conc_format=="log10trns") {
      df$conc = 10^df$conc
    }
    
    # Fit IC50 curves
    min.conc = min(df$conc)
    max.conc = max(df$conc)
    curve.doses = 10^seq(log10(min.conc),log10(max.conc),length.out=100)
    
    cndtns = unique(df$serum)
    
    curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
    curves.df = curves.df[,c(2,1)]
    curves.df$log.dose = log10(curves.df$dose)
    curves.df$pred = NA
    
    IC50s = list(condition = c(T,F),
                 value = c())
    
    for (cndtn in cndtns) {
      
      # Subset on condition
      sub.df = df[df$serum==cndtn,]
      
      # Fit 4-point logistic
      dr4pl.fit = dr4pl(rel.viab~conc,
                        data=sub.df,
                        method.init="logistic")
      # init.parm = dr4pl_theta(0.99,100,-2,0.01),
      # upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
      
      # Predict dose response for curve.doses range
      pred = MeanResponse(coef(dr4pl.fit),curves.df$dose[curves.df$condition==cndtn])
      curves.df$pred[curves.df$condition==cndtn] = pred
      
      # IC50s$value[IC50s$condition==cndtn] = coef(dr4pl.fit)[2] # Save IC50 value
      ## TO DO: update this to absolute IC50
      
      # Calculate the absolute IC50
      C = coef(dr4pl.fit)
      IC50s$value[IC50s$condition==cndtn] = C[2]*((C[4]-C[1])/(0.5-C[1])-1)^(1/C[3]) # Inverse of hill function for y=0.5
      
      
    }
    
    # Make ordering of serum conditions consistent
    df$serum = factor(df$serum,levels=c(T,F))
    curves.df$condition = factor(curves.df$condition,levels=c(T,F))
    
    # Plot colors
    col.srm = "#785DAC"
      col.no.srm = "#999999"
        
      # IC50 values
      IC50.srm = IC50s$value[IC50s$condition==T]
      IC50.no.srm = IC50s$value[IC50s$condition==F]
      
      # Breaks for x axis - remove those closest to IC50 values to avoid overlap
      scale.breaks = 10^(log10(min.conc):log10(max.conc))
      scale.breaks.diff1 = which.min(abs(scale.breaks-IC50.srm))
      scale.breaks.diff2 = which.min(abs(scale.breaks-IC50.no.srm))
      scale.breaks = scale.breaks[-c(scale.breaks.diff1,scale.breaks.diff2)]
      
      
      ggplot()+theme_bw()+
        geom_line(data=curves.df,aes(x=dose,y=pred,color=condition),size=2.5,alpha=0.75)+
        geom_point(data=df,aes(x=conc,y=rel.viab,color=serum),size=5)+
        geom_segment(aes(x=IC50.srm,xend=IC50.srm,y=0.5,yend=0),color=col.srm,size=2,arrow=arrow())+ # Plot IC50 w/o serum
        geom_segment(aes(x=IC50.no.srm,xend=IC50.no.srm,y=0.5,yend=0),color=col.no.srm,size=2,arrow=arrow())+ # Plot IC50 w/ serum
        geom_segment(aes(x=IC50.no.srm,xend=IC50.srm,y=0.25,yend=0.25),color="#087e8b",size=1.5,arrow=arrow(ends="both"))+ # Plot IC50 w/ serum
        scale_x_continuous(breaks=c(scale.breaks,IC50.srm,IC50.no.srm),
                           labels=as.character(round(c(scale.breaks,IC50.srm,IC50.no.srm),1)),
                           trans="log10")+
        xlab(paste0("Drug Concentration (",input$conc_unit,")"))+
        ylab("Relative Viability")+
        scale_color_manual(values=c(col.srm,col.no.srm))+
        theme(
          panel.grid.minor = element_blank(),
          axis.title = element_text(face="bold",size=12),
          axis.text = element_text(face="bold",size=10),
          axis.text.x = element_text(angle = 45,hjust=1,
                                     color=c(rep("black",length(scale.breaks)),col.srm,col.no.srm))
        )+
        guides(color="none")
  })
  
  output$shift_table = renderTable({
    df = dr_df()
    colnames(df) = c("conc","serum","replicate","rel.viab")
    
    # Transform drug concentrations
    if(input$conc_format=="log10trns") {
      df$conc = 10^df$conc
    }
    
    # Fit IC50 curves
    min.conc = min(df$conc)
    max.conc = max(df$conc)
    curve.doses = 10^seq(log10(min.conc),log10(max.conc),length.out=100)
    
    cndtns = unique(df$serum)
    
    curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
    curves.df = curves.df[,c(2,1)]
    curves.df$log.dose = log10(curves.df$dose)
    curves.df$pred = NA
    
    IC50s = list(condition = c(T,F),
                 value = c())
    
    for (cndtn in cndtns) {
      
      # Subset on condition
      sub.df = df[df$serum==cndtn,]
      
      # Fit 4-point logistic
      dr4pl.fit = dr4pl(rel.viab~conc,
                        data=sub.df,
                        method.init="logistic")
      # init.parm = dr4pl_theta(0.99,100,-2,0.01),
      # upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
      
      # Predict dose response for curve.doses range
      pred = MeanResponse(coef(dr4pl.fit),curves.df$dose[curves.df$condition==cndtn])
      curves.df$pred[curves.df$condition==cndtn] = pred
      
      # IC50s$value[IC50s$condition==cndtn] = coef(dr4pl.fit)[2] # Save IC50 value
      ## TO DO: update this to absolute IC50
      
      # Calculate the absolute IC50
      C = coef(dr4pl.fit)
      IC50s$value[IC50s$condition==cndtn] = C[2]*((C[4]-C[1])/(0.5-C[1])-1)^(1/C[3]) # Inverse of hill function for y=0.5
      
      
    }
    
    # Make ordering of serum conditions consistent
    df$serum = factor(df$serum,levels=c(T,F))
    curves.df$condition = factor(curves.df$condition,levels=c(T,F))
    
    # Plot colors
    col.srm = "#785DAC"
      col.no.srm = "#999999"
        
      # IC50 values
      IC50.srm = IC50s$value[IC50s$condition==T]
      IC50.no.srm = IC50s$value[IC50s$condition==F]
      srm.shift = IC50.srm/IC50.no.srm
      
      Cave = as.numeric(input$Cave)
      
      if (!is.na(Cave)) {
        eff.Cave = Cave/srm.shift
        results = data.frame("IC50" = IC50.no.srm,
                             "IC50_shift" = IC50.srm,
                             "serum_shift" = srm.shift,
                             "Cave" = Cave,
                             "eff_Cave" = eff.Cave)
      } else {
        results = data.frame("IC50" = IC50.no.srm,
                             "IC50_shift" = IC50.srm,
                             "serum_shift" = srm.shift)
      }
      
  })
  
  
}
