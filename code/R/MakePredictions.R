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
# stan fits 
best.M1=readRDS(file.path(model.dir, 'M1_4_J_fit.RDS'))
best.M2=readRDS(file.path(model.dir, 'M2_4_Z_fit.RDS'))
# get posteriors 
M1.posteriors=format.posteriors(best.M1)
M2.posteriors=format.posteriors(best.M2)
## index of food items
index.F.M1=read_csv("results/other_tables/prey_index_M1.csv")
index.F.M2=read_csv("results/other_tables/prey_index_M2.csv")
# vectors for z1 and z2
temp.vec=generate.transformation()[[1]]
food.vec=generate.transformation()[[2]]

### predict IGRt ####
igr.prediction=predict.igr(M1.posteriors = M1.posteriors,
            M2.posteriors = M2.posteriors, 
            time.vec = seq(0,220,10),
            experimental.temperature = 19, 
            proportion.libitum = 0.8,
            prey= 'medium_shrimp', # "small_shrimp";"fish";"medium_shrimp";"shrimp_mix"; "amphipods";"starved shrimp";"crayfish"  
            food.vec = food.vec,
            temp.vec = temp.vec) 

ggplot(data=igr.prediction, aes(x=t, y=IGR.mu))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=IGR.lo, ymax=IGR.hi),
              alpha=0.2)+
  theme_minimal()+
  ylab('IGR (%bw*day-1')+
  xlab('Time (DAH)')+
  theme(legend.position = 'bottom')
ggsave(filename='C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph/results/plots/example_plot.jpeg', width = 10, height = 5, units='cm')

### weight
library(pracma)
season.file=data.frame(day.ini=c(1,51,101), 
                       day.fin=c(50,100,150), 
                       temp=c(18,17,16), 
                       proportion.libitum=c(0.8,0.8,0.8),
                       food_type='medium_shrimp') # "small_shrimp";"fish";"medium_shrimp";"shrimp_mix"; "amphipods";"starved shrimp";"crayfish"  
hatch.distr=rnorm(n=2000, mean=0.9, sd=0.1)

w.prediction=predict.integral(M1.posteriors = M1.posteriors,
                 M2.posteriors = M2.posteriors,
                 season.file = season.file,
                 food.vec = food.vec,
                 temp.vec = temp.vec) 

weights=as.data.frame(w.prediction$w)
names(weights)=c('Day1','Day51', 'Day101','Day151')
weights$unit='grams'

ggplot(data = weights)+
  theme_minimal()+
  geom_density(aes(x=Day101))+
  xlab('Grams at 101 DAH')
ggsave(filename='C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph/results/plots/example_plot2.jpeg', width = 10, height = 5, units='cm')

weights.tidy=weights%>%pivot_longer(-unit)
weights.tidy$day=as.numeric(str_remove(weights.tidy$name, 'Day'))

ggplot(data = weights.tidy)+
  geom_boxplot(aes(y=value, x=as.factor(day)))










