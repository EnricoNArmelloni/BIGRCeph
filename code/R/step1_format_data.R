# Description ####
# This script take as input the files st1_experimental details and st_sample_info and applies algorithms described in Appendix_A to estimate missing information. Different methods are applied depending on missing information as described in Table 1. 
rm(list=ls())

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
read.dir='data'
save.plots.dir='results/plots'
setwd(project.dir)

# libraries and functions
library(tidyverse)
library(readxl)
library(mgcv)
source('code/R/supporting_functions.R')

# options
gam.diags=F

# Load tabular data ####
# tabular data are those reported as tables in the publications
experiments.data=read.csv(file.path(read.dir, "st2_sample_info.csv"))
experiments.info=read_delim(file.path(read.dir, "st1_experimental_details.csv"), 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

# exclude non suitable treatments
experiments.data=experiments.data[experiments.data$exp_id %in% experiments.info[experiments.info$included=='y',]$exp_id,]

# paste diet information and polish
diet=read_csv("data/diet_info.csv")
experiments.data=experiments.data%>%left_join(diet, by='exp_id')

## GAM model to handle case III
dat.model=experiments.data[!is.na(experiments.data$igr.sd),]
dat.model$log.igr.sd=log(dat.model$igr.sd)
mod.err=gam(log.igr.sd ~ s(t1,t2,k=10) +  igr.mu + s(nsamp,k=5), data=dat.model)
error.fixed=median(dat.model$igr.sd)

if(gam.diags ==TRUE){
  # diags for gam model
  jpeg(file.path(save.plots.dir, 'gam.jpeg'), width = 20, height = 20, units='cm', res=500)
  par(mfrow = c(2, 2))  # 2x2 plotting layout
  plot(residuals(mod.err) ~ fitted(mod.err), main = "Residuals vs Fitted", 
       xlab = "Fitted Values", ylab = "Residuals")
  abline(h = 0, col = "red", lty = 2)
  hist(residuals(mod.err), breaks = 15, main = "Histogram of Residuals",
       col = "darkgrey", xlab = "Residuals")
  qqnorm(residuals(mod.err), main = "Q-Q Plot")
  qqline(residuals(mod.err), col = "red",  lty = 2)
  plot(dat.model$igr.mu, residuals(mod.err), 
       xlab = expression(IGR[mu]))
  abline(h = 0, col = "red", lty = 2)
  dev.off()
  jpeg(file.path(save.plots.dir, 'gam_effects.jpeg'), width = 20, height = 10, units='cm', res=500)
  par(mfrow = c(1,2))  # 2x2 plotting layout
  plot.gam(mod.err, scheme=2)
  dev.off()
}

## loop over the data and apply estimations when needed ####
x.treatments=unique(experiments.data$exp_id)
x.dat.gr=experiments.data
exp.gr=unique(x.dat.gr$exp_id)
x.dat.gr$case=NA
sample.store=NULL
for(k in 1:nrow(x.dat.gr)){
    cat(k)
    jdat=x.dat.gr[k,]
    case=1
    notes=NA
    
    # detect case. 
    if(is.na(jdat$igr.sd)){
      if(is.na(jdat$w1)|is.na(jdat$w2)){
        case=3
      }else{
        if(is.na(jdat$igr.mu)){
          case=4
        }else{
          case=2
        }
      }
      if(is.na(jdat$igr.mu)&is.na(jdat$igr.sd)&is.na(jdat$errorw2)&is.na(jdat$errorw1)){
        case=5
      }
    }
    
  x.dat.gr[k,]$case=case
  if(case==1){
    # Case 1: No reconstruction needed.
    next
  }
  
  if(case==2){
    # igr.sd is estimated using algorithm EM3.
    w.matrix=data.frame(w.mu.pop = c(x.dat.gr[k,]$w1, x.dat.gr[k,]$w2),
                        w.sd.pop=c(x.dat.gr[k,]$errorw1, x.dat.gr[k,]$errorw2))
    biol.dat=data.frame(t1=x.dat.gr[k,]$t1, t2=x.dat.gr[k,]$t2, igr.mu=x.dat.gr[k,]$igr.mu)
    igr.sd.pred=EM3(w.matrix = w.matrix, biol.value = biol.dat, s.size=50,iter=50,sd.vector=seq(0.1,5,0.2))
    x.dat.gr[k,]$igr.sd=as.numeric(igr.sd.pred[,1])
  }
  
  if(case==3){
    # igr.sd is estimated using GAM model
    x.dat.gr[k,]$igr.sd=exp(predict(mod.err, x.dat.gr[k,]))
  }
  
  if(case==4){
    # igr.mu and igr.sd are estimated from weight-at-age using algorithms EM1 and EM3.
    
    w.matrix=data.frame(w.mu.pop = c(x.dat.gr[k,]$w1, x.dat.gr[k,]$w2),
                        w.sd.pop=c(x.dat.gr[k,]$errorw1, x.dat.gr[k,]$errorw2))
    biol.dat=data.frame(t1=x.dat.gr[k,]$t1, t2=x.dat.gr[k,]$t2, igr.mu=x.dat.gr[k,]$igr.mu)
    est.1=EM1(w.matrix = w.matrix, biol.value = biol.dat)
    biol.dat$igr.mu=est.1$gr.pred
    est.2=EM3(w.matrix = w.matrix, biol.value = biol.dat, iter=50, s.size=100, sd.vector = seq(0.1,5,0.1))
    
    x.dat.gr[k,]$igr.mu=est.1$gr.pred
    x.dat.gr[k,]$igr.sd=as.numeric(est.2[,1])
  }
  
  if(case==5){
    # igr.mu is estimated from weight-at-age using algorithm EM1, igr.sd is estimated using the GAM.
    
    w.matrix=data.frame(w.mu.pop = c(x.dat.gr[k,]$w1, x.dat.gr[k,]$w2),
                        w.sd.pop=c(0.1,0.1))
    biol.dat=data.frame(t1=x.dat.gr[k,]$t1, t2=x.dat.gr[k,]$t2, igr.mu=x.dat.gr[k,]$igr.mu)
    est.1=EM1(w.matrix = w.matrix, biol.value = biol.dat)
    x.dat.gr[k,]$igr.mu=est.1$gr.pred
    x.dat.gr[k,]$igr.sd=exp(predict(mod.err, x.dat.gr[k,]))
  }
}

x.dat.gr=x.dat.gr[!is.na(x.dat.gr$igr.mu),]

# inspect the number of data per case and aggregation
x.dat.gr%>%distinct(publication_id)
x.dat.gr%>%distinct(experiment)
x.dat.gr%>%distinct(exp_id)

x.dat.gr%>%
  dplyr::filter(t<=120)%>%
  #distinct(exp_id, notes)%>%
  dplyr::group_by(case)%>%
  tally()

# number of treatmens
x.dat.gr%>%
  dplyr::filter(t<=120)%>%
  distinct(exp_id, case)%>%
  dplyr::group_by(case)%>%
  tally()
# number of publications
x.dat.gr%>%
  dplyr::filter(t<=120)%>%
  dplyr::mutate(publ=substr(experiment, 1, nchar(experiment)-2))%>%
  distinct(publ, case)%>%
  dplyr::group_by(case)%>%
  tally()

x.dat.gr%>%
  dplyr::filter(t>120)%>%
  #distinct(exp_id, notes)%>%
  dplyr::group_by(case)%>%
  tally()
# number of treatmens
x.dat.gr%>%
  dplyr::filter(t>120)%>%
  distinct(exp_id, case)%>%
  dplyr::group_by(case)%>%
  tally()
# number of publications
x.dat.gr%>%
  dplyr::filter(t>120)%>%
  dplyr::mutate(publ=substr(experiment, 1, nchar(experiment)-2))%>%
  distinct(publ, case)%>%
  dplyr::group_by(case)%>%
  tally()

# save
write.csv(x.dat.gr, file.path(read.dir, 'input_data.csv'), row.names = F)

