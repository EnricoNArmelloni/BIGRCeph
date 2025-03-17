# Description ####
# this script demonstrate the use of MCMC draws to predict IGRt 
rm(list=ls())

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
read.dir='results'
model.dir='results/model_fits'
save.dir='results'
setwd(project.dir)

# libraries and functions
library(tidyverse)
source('code/R/supporting_functions.R')

### load data ####
### get stan fits 
best.M1=readRDS(file.path(model.dir, 'M1_4_J_fit.RDS'))
best.M2=readRDS(file.path(model.dir, 'M2_4_Z_fit.RDS'))

### get input data
M1.data=read_csv("data/M1_model_data.csv")
M2.data=read_csv("data/M2_model_data.csv")

### get posteriors 
M1.posteriors=format.posteriors(best.M1)
M2.posteriors=format.posteriors(best.M2)

## match fits index with food items
index.F.M1=distinct(M1.data, Index.food, main_prey)
index.F.M2=distinct(M2.data, Index.food, main_prey)


### predict IGRt ####
# vectors for z1 and z2
temp.vec=data.frame(temp=seq(12.5,27, length.out=1000),
                    z2=seq(0, 1, length.out=1000))
temp.vec$temp=round(temp.vec$temp, 2)

food.vec=data.frame(Rcomma=seq(0,1.66, length.out=166),
                    z1=seq(0, 1, length.out=166))
food.vec$Rcomma=round(food.vec$Rcomma, 2)

# prediction 
# define settings
time.vec=seq(0,120,5)
experimental.temperature=18
proportion.libitum=0.8
prey='medium_shrimp' # "small_shrimp";"fish";"medium_shrimp";"shrimp_mix"; "amphipods";"starved shrimp";"crayfish"  

# predictions for M1 and M2
m1.store=matrix(nrow = nrow(M1.posteriors), ncol=length(time.vec))
m2.store=matrix(nrow = nrow(M2.posteriors), ncol=length(time.vec))
for(i in 1:length(time.vec)){  
  m1.pred=data.frame(pred=M1(a=M1.posteriors$a,
                               b1=M1.posteriors$bF,
                               b2=M1.posteriors$bT1,
                               b3=M1.posteriors$bT2,
                               sigma=M1.posteriors$`sigma[28]`,
                               mu=M1.posteriors$mu,
                               t=time.vec[i],
                               tem=temp.vec[temp.vec$temp==experimental.temperature,]$z2,
                               Fi=M1.posteriors[,paste0('F[',index.F.M1[index.F.M1$main_prey==prey,]$Index.food,']')],
                               R=food.vec[food.vec$Rcomma==proportion.libitum,]$z1))
  m1.pred$pred=as.numeric(mapply(function(mu, sigma) rnorm(n = 1, mean = mu, sd = sigma), m1.pred$pred, 
                                   apply(M1.posteriors[,grep('obs_erro', names(M1.posteriors))],1,mean)))
  m1.store[,i]=m1.pred$pred
    
  m2.pred=data.frame(pred=M2(A=M2.posteriors$A,
                               b1=M2.posteriors$bF,
                               b4a=M2.posteriors$bTA,
                               b2=M2.posteriors$bTB1,
                               b3=M2.posteriors$bTB2,
                               B=M2.posteriors$`B[39]`,
                               C=M2.posteriors$C,
                               t=time.vec[i],
                               tem=temp.vec[temp.vec$temp==experimental.temperature,]$z2,
                               Fi=M2.posteriors[,paste0('F[',index.F.M2[index.F.M2$main_prey==prey,]$Index.food,']')],
                               R=food.vec[food.vec$Rcomma==proportion.libitum,]$z1))
  m2.pred$pred=as.numeric(mapply(function(mu, sigma) rnorm(n = 1, mean = mu, sd = sigma), m2.pred$pred, 
                                   apply(M2.posteriors[,grep('obs_erro', names(M2.posteriors))],1,mean)))
  m2.store[,i]=m2.pred$pred
    
  }
m1.store=as.data.frame(m1.store)
names(m1.store)=paste0('t', time.vec)
m2.store=as.data.frame(m2.store)
names(m2.store)=paste0('t', time.vec)

## apply sigmoid trnasition
m.comb=as.data.frame(apply_sigmoid_row(m1.store,m2.store,time.vec))

# format for plotting
Mcomb.resp=m.comb%>%
  dplyr::mutate(Temperature=experimental.temperature)%>%
  pivot_longer(-c(Temperature), names_to = 't', values_to = 'IGR')%>%
  dplyr::mutate(t=as.numeric(str_remove(t,'t')))%>%
  dplyr::group_by(t,Temperature)%>%
  dplyr::summarise(IGR.mu=mean(IGR, probs=c(0.5)),
                   IGR.lo=quantile(IGR, probs=c(0.1)),
                   IGR.hi=quantile(IGR, probs=c(0.9)))


ggplot(data=Mcomb.resp, aes(x=t, y=IGR.mu))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=IGR.lo, ymax=IGR.hi),
              alpha=0.2)+
  theme_minimal()+
  ylab('IGR (%bw*day-1')+
  xlab('Time (DAH)')+
  theme(legend.position = 'bottom')

ggsave(filename='C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph/results/plots/example_plot.jpeg', width = 10, height = 5, units='cm')
