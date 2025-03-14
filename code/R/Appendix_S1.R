rm(list=ls())
# Description ####
# Description ####
# We aim to test different strategies to estimate mean and uncertainty of growth rate starting from aggregated values. We are interested in estimating both the effect of the sample size and of .... We test three strategies: population statistic, resampling and sensitivity

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
read.dir='data'
save.dir='results/diags'
setwd(project.dir)


# libraries and functions ####
library(tidyverse)
library(readxl)
library(ggExtra)
library(mgcv)
library(ggridges)
library(doParallel)
source('code/R/supporting_functions.R')

# options
theme_set(theme_bw())


## load real data ####
experiments.data=read.csv(file.path(read.dir, "st2_sample_info.csv"))
experiments.info=read_delim(file.path(read.dir, "st1_experimental_details.csv"), 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

# exclude data
experiments.data=experiments.data[experiments.data$exp_id %in% experiments.info[experiments.info$included=='y',]$exp_id,]


## derive quantities to start the population ####
# hatchlings
hatch.dat=experiments.data[experiments.data$t1<=1 & !is.na(experiments.data$w1),]%>%
  distinct(publication_id, t1,w1,errorw1 )%>%
  dplyr::group_by(publication_id)%>%
  dplyr::summarise(w1=mean(w1), errorw1=mean(errorw1))

## define vector of mean GR
gr.vector=data.frame(igr.mu=c(4,9.5,8.4,6.5,4.6,3.5,4.3,3,2,1.8,1.6,2,1.4,1.1,1),
                     t1=NA,
                     t2=seq(11,151,10))
gr.vector$t1=lag(gr.vector$t2)
gr.vector$t1[1]=1
gr.vector$t=(gr.vector$t1+gr.vector$t2)/2

p.mu.gr=ggplot(data=experiments.data[experiments.data$t2<=150,], aes(x=((t)), y=igr.mu))+
  geom_point()+
  geom_line(data=gr.vector, color='red')+
  geom_point(data=gr.vector, color='red')+
  xlab('Time (Days)')+
  ylab(expression(mu[IGR][t]));p.mu.gr

## gr sd
plot1.dat=experiments.data
plot1.dat$log.sd=log(experiments.data$igr.sd)
plot1.dat$period=ifelse(plot1.dat$t<=21, 'Period1',
                        ifelse(plot1.dat$t>61, 'Period3','Period2'))

p1=ggplot(data=plot1.dat[plot1.dat$t2<150,],aes(x=t, y=(igr.sd), color=period))+
  geom_point()+
  geom_vline(xintercept = 21)+
  #geom_vline(xintercept = 31)+
  geom_vline(xintercept = 61)+
  theme(legend.position = 'bottom')+ylab(expression(sigma[IGR][t]))+
  xlab('Time (Days)')+
  scale_color_manual(values=c("#31688EFF","#440154FF", "#6DCD59FF" ))

p.sd.gr=ggMarginal(p1, groupColour = TRUE, groupFill = TRUE, margins = 'y')
fig.2=ggpubr::ggarrange(p.mu.gr,p.sd.gr, ncol=1,  labels = c("a)", "b)"))
ggsave('results/plots/Appendix_S1/observed_estimated_hatch_w.png', width = 15, height = 10, units='cm', dpi=500)

experiments.data=experiments.data[!is.na(experiments.data$igr.sd),]
set.seed(46)
gr.distribution.1=lnorm.transformation(mean(experiments.data[experiments.data$t<=21,]$igr.sd),
                                       sd(experiments.data[experiments.data$t<=21,]$igr.sd), out='distr', npop=500);gr.distribution.1
set.seed(46)
lnorm.transformation(mean(experiments.data[experiments.data$t<=21,]$igr.sd),sd(experiments.data[experiments.data$t<=21,]$igr.sd),npop=500)

set.seed(46)
gr.distribution.2=lnorm.transformation(mean(experiments.data[experiments.data$t%in%22:60,]$igr.sd),
                                       sd(experiments.data[experiments.data$t%in%22:60,]$igr.sd), out='distr', npop=500)
set.seed(46)
lnorm.transformation(mean(experiments.data[experiments.data$t%in%22:60,]$igr.sd),sd(experiments.data[experiments.data$t%in%22:60,]$igr.sd),  npop=500)

set.seed(46)
gr.distribution.3=lnorm.transformation(mean(experiments.data[experiments.data$t>60,]$igr.sd),
                                       sd(experiments.data[experiments.data$t>60,]$igr.sd), out='distr', npop=500)

### experiment
npop=1000
nsamp=30
nsims=100

for(j.population in 1:100){
  
  cat('starting', j.population)
  
  set.seed(j.population)
  # population parameters
  sd.vec=c(sample(gr.distribution.1,2),
         sample(gr.distribution.2,4),
         sample(gr.distribution.3,9))
  biol.pars=gr.vector
  biol.pars$igr.sd=sd.vec

  # generate population
  pop.matrix=matrix(nrow=npop, ncol=nrow(biol.pars)+1)
  pop.matrix[,1]=lnorm.transformation(mean(hatch.dat$w1), mean(hatch.dat$errorw1), out='distr', npop = npop)
  gr.matrix=matrix(nrow=npop, ncol=nrow(biol.pars))
  for(i in 1:nrow(biol.pars)){
    gr.i=rnorm(n=npop, mean=biol.pars[i,]$igr.mu, sd=biol.pars[i,]$igr.sd)
    pop.matrix[,i+1]=un.gr.formula(w1=pop.matrix[,i], gr=gr.i, t1=biol.pars[i,]$t1,t2=biol.pars[i,]$t2)
    gr.matrix[,i]=gr.i
  }

  ## test estimation strategies on the entire population ####
  w.means=data.frame(w.mu.pop=apply(pop.matrix,2,mean),
                     w.sd.pop=apply(pop.matrix,2,sd))
  EM1.predictions.pop=EM1(w.matrix = w.means, biol.value = biol.pars)
  EM2.prediction.pop=EM2(w.matrix = w.means, biol.value = biol.pars, s.size=500, iter=50)
  biol.pars.pred=biol.pars
  biol.pars.pred$igr.mu=EM2.prediction.pop$mu.mu
  EM3.predictions.pop=EM3(w.matrix = w.means, biol.value = biol.pars.pred, s.size=500, iter=50,sd.vector=seq(0.1,5,0.2))

  outputs.pop=data.frame(t=biol.pars$t,
             igr.mu.obs=biol.pars$igr.mu,
             igr.sd.obs=round(biol.pars$igr.sd, digits=3),
             igr.mu.pred.EM1=round(EM1.predictions.pop$gr.pred, digits=2),
             igr.sd.pred.EM1=round(EM1.predictions.pop$igr.sd.pred, digits=3),
             igr.mu.pred.EM2=round(EM2.prediction.pop$mu.mu, digits=2),
             igr.mu.pred.EM2.dev=round(EM2.prediction.pop$mu.sd, digits=2),
             igr.sd.pred.EM2=round(EM2.prediction.pop$sd.mu, digits=3),
             igr.sd.pred.EM2.dev=round(EM2.prediction.pop$sd.sd.s, digits=3),
             igr.sd.pred.EM3=round(EM3.predictions.pop[,1], digits=3))
  
  ## test estimation strategies on a sample from the population ####
  niter=30
  xx=2
  cl <- parallel::makeCluster(10)
  doParallel::registerDoParallel(cl)
  time.ini=Sys.time()
  stores=foreach(xx = 1:15,
          .packages = c("tidyverse"))%dopar%{
    store.sample.results=list()
    for(xx in 1:2){
    cat(xx)
      
    # sample from population and eliminate sampled animal
    sample.matrix=matrix(nrow=nsamp, ncol=nrow(biol.pars)+1) # store here the weights of the sampled individuals
    sample.gr=matrix(nrow=nsamp, ncol=nrow(biol.pars)) # store here the real GR associated to the sample
    pop.l=1:npop
    for(i in 1:ncol(pop.matrix)){
      sampled=sample(pop.l,nsamp,replace = F) # sample w
      sample.matrix[,i]=pop.matrix[sampled,i]
      if(i<ncol(pop.matrix)){
        sample.gr[,i]=gr.matrix[sampled,i]   # sample ge
      }
      pop.l=pop.l[pop.l%in%sampled==F] # remove sampled individuals
    }
    
    # compute aggregated values
    w.means=data.frame(w.mu.pop=apply(sample.matrix,2,mean),
                       w.sd.pop=apply(sample.matrix,2,sd))
    gr.sample.obs=data.frame(igr.mu.samp=apply(sample.gr,2,mean),
                             igr.sd.samp=apply(sample.gr,2,sd))
    
  
  
    ## test estimation strategies on the sample ####
    EM1.predictions.samp=EM1(w.matrix = w.means, biol.value = biol.pars)
    EM2.prediction.samp=EM2(w.matrix = w.means, biol.value = biol.pars, s.size=30, iter=niter)
    biol.pars.pred=biol.pars
    biol.pars.pred$igr.mu=EM2.prediction.samp$mu.mu
    EM3.predictions.samp=EM3(w.matrix = w.means, biol.value = biol.pars.pred, s.size=30, iter=niter,sd.vector=seq(0.1,5,0.2))
    
    outputs.samp=data.frame(t=biol.pars$t,
                           igr.mu.obs=biol.pars$igr.mu,
                           igr.mu.sample=gr.sample.obs$igr.mu.samp,
                           igr.sd.obs=round(biol.pars$igr.sd, digits=3),
                           igr.sd.sample=round(gr.sample.obs$igr.sd.samp, digits=3),
                           igr.mu.pred.EM1=round(EM1.predictions.samp$gr.pred, digits=2),
                           igr.sd.pred.EM1=round(EM1.predictions.samp$igr.sd.pred, digits=3),
                           igr.mu.pred.EM2=round(EM2.prediction.samp$mu.mu, digits=2),
                           igr.mu.pred.EM2.dev=round(EM2.prediction.samp$mu.sd, digits=2),
                           igr.sd.pred.EM2=round(EM2.prediction.samp$sd.mu, digits=3),
                           igr.sd.pred.EM2.dev=round(EM2.prediction.samp$sd.sd.s, digits=3),
                           igr.sd.pred.EM3=round(EM3.predictions.samp[,1], digits=3))
    # store information
    store.sample.results[[xx]]=outputs.samp
    }
    
    # END of loop: rearrange results and return it
    store.sample.results=plyr::ldply(store.sample.results)
    return(store.sample.results)
  }
  
  time.fin=Sys.time()
  print(time.fin-time.ini)
  stopCluster(cl)
  
  # disaggregate results
  store.sample.results.tidy=plyr::ldply(stores)
  write.csv(outputs.pop,paste0('results/Appendix_S1/estimates_pop_rep_',j.population ,'.csv'),row.names = F)
  write.csv(store.sample.results.tidy,paste0('results/Appendix_S1/estimates_sample_rep_',j.population ,'.csv'),row.names = F)
  
}

## performance metrics ####
x.res=list.files(path='results/Appendix_S1', pattern='.csv')
x.res=x.res[grep('sample', x.res)]

store.sig=NULL
for(i in 1:length(x.res)){
  i.res=read.csv(file.path('results/Appendix_S1', x.res[i]))
  i.res=data.frame(t=unique(i.res$t),
  Y1.mu=as.numeric(testing.predictions(xdat=i.res, xobs='igr.mu.sample', xvar='igr.mu.pred.EM1', alpha=0.01)),
  Y2.mu=as.numeric(testing.predictions(xdat=i.res, xobs='igr.mu.sample', xvar='igr.mu.pred.EM2', alpha=0.01)),
  Y1.sd=as.numeric(testing.predictions(xdat=i.res, xobs='igr.sd.sample', xvar='igr.sd.pred.EM1', alpha=0.01)),
  Y2.sd=as.numeric(testing.predictions(xdat=i.res, xobs='igr.sd.sample', xvar='igr.sd.pred.EM2', alpha=0.01)),
  Y3.sd=as.numeric(testing.predictions(xdat=i.res, xobs='igr.sd.sample', xvar='igr.sd.pred.EM3', alpha=0.01)))
  i.res$iter=i
  store.sig=rbind(store.sig, i.res)
}


table.GRsd.t=store.sig%>%
  dplyr::group_by(t)%>%
  dplyr::summarise(e1=round(1-(mean(Y1.sd)), digits=3),
                   e2=round(1-(mean(Y2.sd)),digits=3),
                   e3=round(1-(mean(Y3.sd)), digits=3))
table.GRsd.average=store.sig%>%
  dplyr::summarise(e1=round(1-(mean(Y1.sd)), digits=3),
                   e2=round(1-(mean(Y2.sd)),digits=3),
                   e3=round(1-(mean(Y3.sd)), digits=3))
table.GRmu.t=store.sig%>%
  dplyr::group_by(t)%>%
  dplyr::summarise(e1=round(1-(mean(Y1.mu)), digits=3),
                   e2=round(1-(mean(Y2.mu)),digits=3))
table.GRmu.average=store.sig%>%
  dplyr::summarise(e1=round(1-(mean(Y1.mu)), digits=3),
                   e2=round(1-(mean(Y2.mu)),digits=3))

write.csv(table.GRmu.t, 'results/sim_table_1_t.csv', row.names = F)
write.csv(table.GRmu.average, 'results/sim_table_1_avg.csv', row.names = F)
write.csv(table.GRsd.t, 'results/sim_table_2_t.csv', row.names = F)
write.csv(table.GRsd.average, 'results/sim_table_2_avg.csv', row.names = F)


## plots ####
# population
pop.res=read_csv("results/Appendix_S1/estimates_pop_rep_10.csv")
pop.res.mu=pop.res[,c('t','igr.mu.obs','igr.mu.pred.EM1','igr.mu.pred.EM2','igr.mu.pred.EM2.dev')]
names(pop.res.mu)[3:4]=c('EM1','EM2')
pop.res.mu=pop.res.mu%>%pivot_longer(-c('t','igr.mu.obs' ,'igr.mu.pred.EM2.dev'), names_to = 'method', values_to = 'gr.pred')

p0=ggplot(data=pop.res.mu)+
  geom_point(aes(x=t,y=gr.pred, color=method))+
  theme(legend.position = 'bottom')+
  ylab(expression(mu[IGR][t]))+
  xlab('Time (DAH)')+
  labs(color='Method')+
  scale_color_manual(values=c("#31688EFF","#440154FF", "#6DCD59FF" ), labels=c('EM2','EM1'))+
  geom_errorbar(aes(x=t,ymin=gr.pred-igr.mu.pred.EM2.dev, ymax=gr.pred+igr.mu.pred.EM2.dev, color=method))+
  geom_point(aes(x=t, y=igr.mu.obs))

pop.res.sd=pop.res[,c('t','igr.sd.obs','igr.sd.pred.EM1','igr.sd.pred.EM2','igr.sd.pred.EM2.dev','igr.sd.pred.EM3')]
names(pop.res.sd)[c(3,4,6)]=c('EM1','EM2','EM3')
#pop.res.sd=pop.res.sd%>%pivot_longer(-c('t','igr.sd.obs' ,'igr.sd.pred.EM2.dev'), names_to = 'method', values_to = 'gr.pred')

p1=ggplot(data=pop.res.sd)+
  ylab('GR sd')+
  xlab('Time (DAH)')+
  geom_point(aes(x=t+0.5, y=EM2), color='#31688EFF')+
  geom_line(aes(x=t+0.5, y=EM2), color='#31688EFF')+
  geom_errorbar(aes(x=t+0.5, ymin=EM2-igr.sd.pred.EM2.dev,
                    ymax=EM2+igr.sd.pred.EM2.dev), color='#31688EFF')+
  geom_point(aes(x=t, y=igr.sd.obs),shape=18)+
  ylab(expression(sigma[IGR][t]))+
  ggtitle('EM2');p1

p2=ggplot(data=pop.res.sd)+
  ylab(expression(sigma[IGR][t]))+
  xlab('Time (DAH)')+
  geom_line(aes(x=t, y=EM3), color='#6DCD59FF')+
  geom_point(aes(x=t, y=EM3), color='#6DCD59FF')+
  geom_point(aes(x=t, y=igr.sd.obs),shape=18)+
  ggtitle('EM3');p2

p3=ggplot(data=pop.res.sd)+
  ylab(expression(sigma[IGR][t]))+
  xlab('Time (DAH)')+
  geom_line(aes(x=t, y=EM1), color='#440154FF')+
  geom_point(aes(x=t, y=EM1), color='#440154FF')+
  geom_point(aes(x=t, y=igr.sd.obs),shape=18)+
  ggtitle('EM1');p3

p4=ggpubr::ggarrange(p0,p3,p1,p2, labels=c('a)','b)','c)','d)'))  
ggsave(plot=p4, paste0('results/plots/Appendix_S1/simulation_population_rep_',j.population ,'.jpeg'), width = 25, height = 20, units='cm')  


# sample
samp.res=read_csv("results/Appendix_S1/estimates_sample_rep_10.csv")
samp.res.mu=samp.res%>%dplyr::group_by(t)%>%
  dplyr::summarise(EM1=mean(igr.mu.pred.EM1),
                   EM1.sd=sd(igr.mu.pred.EM1),
                   EM2=mean(igr.mu.pred.EM2),
                   EM2.sd=sd(igr.mu.pred.EM2),
                   gr.obs=mean(igr.mu.obs))

p0.dat.p0=samp.res.mu[,c('t','gr.obs', 'EM1', 'EM1.sd')]%>%
  dplyr::mutate(method='EM1')
names(p0.dat.p0)[3:4]=c('pred.mu','pred.sd')
p0.dat.p1=samp.res.mu[,c('t','gr.obs', 'EM2', 'EM2.sd')]%>%
  dplyr::mutate(method='EM2')
names(p0.dat.p1)[3:4]=c('pred.mu','pred.sd')
p0.dat=rbind(p0.dat.p1,p0.dat.p0)

p0=ggplot(data=p0.dat)+
  
  geom_point(aes(x=t-0.5, y=pred.mu, color=method))+
  geom_errorbar(aes(x=t-0.5, ymin=pred.mu-pred.sd,
                    ymax=pred.mu+pred.sd, color=method))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  ylab(expression(mu[IGR][t]))+
  xlab('Time (DAH)')+
  labs(color='Method')+
  scale_color_manual(values=c("#31688EFF", "black","#440154FF" ), labels=c('EM2','EM1', 'EM1'))+
  #geom_point(aes(x=t,y=gr.obs))+
  geom_point(aes(x=t, y=gr.obs),shape=18)


samp.res.sd=samp.res%>%dplyr::group_by(t)%>%
  dplyr::summarise(EM1=mean(igr.sd.pred.EM1),
                   EM1.sd=sd(igr.sd.pred.EM1),
                   EM2=mean(igr.sd.pred.EM2),
                   EM2.sd=sd(igr.sd.pred.EM2),
                   EM3=mean(igr.sd.pred.EM3),
                   EM3.sd=sd(igr.sd.pred.EM3),
                   real.sd=mean(igr.sd.obs),
                   sample.mu=mean(igr.sd.sample),
                   sample.sd=sd(igr.sd.sample))

plot.sample.sd=function(xdat, x.title, x.color, x.var){
xdat=cbind(xdat[, c("t", "real.sd",'sample.sd','sample.mu')],
           xdat[,grep(x.var, names(xdat))])
names(xdat)[5:6]=c('pred.mu','pred.sd')
p=ggplot(data=xdat)+
  ylab(expression(sigma[IGR][t]))+
  xlab('Time (DAH)')+
  
  geom_point(aes(x=t+0.5, y=pred.mu), color=x.color)+
  geom_errorbar(aes(x=t+0.5, ymin=pred.mu-pred.sd,
                    ymax=pred.mu+pred.sd), color=x.color)+
  geom_point(aes(x=t, y=sample.mu),shape=18)+
  geom_point(aes(x=t, y=real.sd),shape=18, color='gold')+
  geom_errorbar(aes(x=t, ymin=sample.mu-sample.sd,
                    ymax=sample.mu+sample.sd), color='black')+
  ggtitle(x.title)
return(p)
}

p1=plot.sample.sd(xdat=samp.res.sd, x.title='EM2', x.color='#31688EFF',x.var='EM2')
p2=plot.sample.sd(xdat=samp.res.sd, x.title='EM3', x.color='#6DCD59FF',x.var='EM3')
p3=plot.sample.sd(xdat=samp.res.sd, x.title='EM1', x.color='#440154FF',x.var='EM1')
p4=ggpubr::ggarrange(p0,p3, p1,p2, labels=c('a)','b)','c)','d)'));p4

ggsave(plot=p4, paste0('results/plots/Appendix_S1/simulation_sample_rep_',j.population ,'.jpeg'), 
       width = 25, height = 20, units='cm')  

























