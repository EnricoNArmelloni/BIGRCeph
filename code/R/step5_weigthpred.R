# Description ####
# this script demonstrate the use of MCMC draws to predict IGRt 
rm(list=ls())
library(pracma)
library(stats4)
library(tidyverse)
library(bayesplot)
library(parallel)
library(foreach)
library(doParallel)

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
read.dir='results'
model.dir='results/model_fits'
save.dir='results'
setwd(project.dir)

# libraries and functions
source('code/R/supporting_functions.R')
color_scheme_set("viridis")
theme_set(theme_bw())

# parameters
nsamp=2000

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

# wild individual data
obs.data=read_csv("data/wild_obs.csv")
xyrs=unique(obs.data$year)
obs.long=NULL
set.seed(46)
for(i in 1:length(xyrs)){
  i.obs=obs.data[obs.data$year==xyrs[i],]
  i.w=base::sample(i.obs$wei_g,nsamp,replace=T, prob=i.obs$Freq)
  i.res=data.frame(year=xyrs[i], weig_g=i.w)
  obs.long=rbind(obs.long, i.res)
}

obs.long%>%
  ggplot(aes(x=weig_g))+
  geom_histogram()+
  facet_wrap(~year)

# hatchlings from Sepe et al 2025
hatch.adr.mu=c(0.08,0.09,0.09,0.08,0.08,0.08)
hatch.adr.sd=c(0.01,0.02,0.02,0.01,0.01,0.01)
egg.w.mu=round(mean(hatch.adr.mu), digits=3)
egg.w.sd=round(sqrt(mean(hatch.adr.sd^2)), digits=3) # pooled sd
set.seed(46)
hatch.distr=lnorm.transformation(mu.w=egg.w.mu, sd.w=egg.w.sd, out='distr', npop=nsamp)


# bottom temperature data ####
obs.temp=read.csv('data/adr_seatemperature.csv')

ggplot(obs.temp, aes(x=month, y=mu.t, color=factor(year)))+
  geom_line()+
  geom_point()+
  facet_wrap(~area)

obs.temp=rbind(obs.temp[obs.temp$area=='coastal'&obs.temp$month%in% 5:9,],
      obs.temp[obs.temp$area=='offshore'&obs.temp$month%in% 10:11,])%>%
  arrange(year,month)

temp.vec=data.frame(temp=seq(12.50,27, 0.01))
temp.vec$z2=seq(0, 1, length.out=nrow(temp.vec))

obs.temp$z2=NA
for(i in 1:nrow(obs.temp)){
  obs.temp[i,]$z2=temp.vec[abs(temp.vec$temp - obs.temp[i,]$mu.t) == min(abs(temp.vec$temp-obs.temp[i,]$mu.t)), ]$z2
}
obs.temp$day.ini=lubridate::yday(as.Date(paste(obs.temp$year, obs.temp$month, '01', sep='-')))-lubridate::yday(as.Date(paste(obs.temp$year, '05', '01', sep='-')))
obs.temp$day.fin=lubridate::yday(as.Date(paste(obs.temp$year, obs.temp$month+1, '01', sep='-')))-lubridate::yday(as.Date(paste(obs.temp$year, '05', '01', sep='-')))

## find best R
R.vec=seq(0.66,0.73,0.005)
yr.res=list()
cl <- makeCluster(length(xyrs))        # Use 4 cores
registerDoParallel(cl)
res.test=foreach(x.yr = 1:length(xyrs) ,.packages = c("pracma"))%dopar%{
  Obs=obs.long[obs.long$year==xyrs[x.yr],]$weig_g
  store.res=NULL
  #x.T=obs.temp.seas[obs.temp.seas$year==xyrs[x.yr],]
  x.T=obs.temp[obs.temp$year==xyrs[x.yr],]
  store.ll=numeric(length = length(R.vec))
  for(x.val in 1:length(R.vec)){
    cat(x.val)
    x.R=R.vec[x.val]
    pop.season=matrix(nrow=nsamp, ncol=nrow(x.T)+1)
    pop.season[,1]=hatch.distr
    gr.season=matrix(nrow=nsamp, ncol=nrow(x.T))
    for(t in 1:nrow(x.T)){
      gr.est=quadv(mean.integral, a=x.T[t,]$day.ini, b=x.T[t,]$day.fin, 
                   aa=M1.posteriors$a,
                   bF1=M1.posteriors$bF,
                   bT1=M1.posteriors$bT1,
                   bT2=M1.posteriors$bT2,
                   mu=M1.posteriors$mu,
                   sigma=M1.posteriors$`sigma[28]`,
                   R=x.R,
                   Fi1=M1.posteriors$`F[4]`,
                   Fi2=M2.posteriors$`F[1]`,
                   tem=x.T[t,]$z2,
                   A= M2.posteriors$A,
                   bTA= M2.posteriors$bTA,
                   bF2= M2.posteriors$bF,
                   B=M2.posteriors$`B[39]`,
                   bTB1= M2.posteriors$bTB1,
                   bTB2= M2.posteriors$bTB2,
                   sigma_obs1= M1.posteriors$sigma_obs,
                   sigma_obs2= M2.posteriors$sigma_obs,
                   C= M2.posteriors$C)
      #gr.vec=gr.est$Q/(day.fin-day.ini)
      gr.season[,t]=gr.est$Q/(x.T[t,]$day.fin-x.T[t,]$day.ini)
      pop.season[,t+1]=un.gr.formula(w1=pop.season[,t], 
                                     gr=gr.season[,t], 
                                     t1=x.T[t,]$day.ini, 
                                     t2=x.T[t,]$day.fin) 
      }
    Pred=pop.season[,t+1]
    NegLogLik <- -1*sum(dnorm(log(Obs),log(Pred),log=T)) 
    store.ll[x.val]=NegLogLik
  }
   x.res=data.frame(nll=store.ll, r.proxy=R.vec, year=xyrs[x.yr])
   yr.res[[x.yr]]=x.res
}
stopCluster(cl)

plyr::ldply(res.test)%>%
  dplyr::group_by(year, r.proxy)%>%
  summarise(nll.mu=mean(nll), nll.sd=sd(nll))%>%
  ggplot(aes(x=r.proxy, y=nll.mu))+
  geom_point()+
  #geom_errorbar(aes(ymin=nll.mu-nll.sd, ymax=nll.mu+nll.sd))+
  facet_wrap(~year)

best.R=plyr::ldply(res.test)%>%
  dplyr::group_by(year, r.proxy)%>%
  summarise(nll.mu=mean(nll), nll.sd=sd(nll))%>%
  dplyr::group_by(year)%>%
  slice_min(nll.mu)
best.R=mean(best.R$r.proxy) # 0.684375

# implementation
yr.pop=list()
cl <- makeCluster(length(xyrs))        # Use 4 cores
registerDoParallel(cl)
pred.pop=foreach(x.yr = 1:length(xyrs) ,.packages = c("pracma"))%dopar%{
  x.T=obs.temp[obs.temp$year==xyrs[x.yr],]
  x.R=best.R
  pop.season=matrix(nrow=nsamp, ncol=nrow(x.T)+1)
  pop.season[,1]=hatch.distr
  gr.season=matrix(nrow=nsamp, ncol=nrow(x.T))
  
  for(t in 1:nrow(x.T)){
    gr.est=quadv(mean.integral, a=x.T[t,]$day.ini, b=x.T[t,]$day.fin, 
                 aa=M1.posteriors$a,
                 bF1=M1.posteriors$bF,
                 bT1=M1.posteriors$bT1,
                 bT2=M1.posteriors$bT2,
                 mu=M1.posteriors$mu,
                 sigma=M1.posteriors$`sigma[28]`,
                 R=x.R,
                 Fi1=M1.posteriors$`F[4]`,
                 Fi2=M2.posteriors$`F[1]`,
                 tem=x.T[t,]$z2,
                 A= M2.posteriors$A,
                 bTA= M2.posteriors$bTA,
                 bF2= M2.posteriors$bF,
                 B=M2.posteriors$`B[39]`,
                 bTB1= M2.posteriors$bTB1,
                 bTB2= M2.posteriors$bTB2,
                 sigma_obs1= M1.posteriors$sigma_obs,
                 sigma_obs2= M2.posteriors$sigma_obs,
                 C= M2.posteriors$C)
    #gr.vec=gr.est$Q/(day.fin-day.ini)
    gr.season[,t]=gr.est$Q/(x.T[t,]$day.fin-x.T[t,]$day.ini)
    pop.season[,t+1]=un.gr.formula(w1=pop.season[,t], 
                                   gr=gr.season[,t], 
                                   t1=x.T[t,]$day.ini, 
                                   t2=x.T[t,]$day.fin) 
    
  }
  yr.pop[[x.yr]]=pop.season
}
stopCluster(cl)

names(pred.pop)=xyrs
pred.weights=plyr::ldply(pred.pop)
names(pred.weights)=c('year', paste0('M',0:7))


### plotting
w.comp=rbind(data.frame(w=obs.long$weig_g, year=obs.long$year, source='Observation'),
             data.frame(w=pred.weights$M7, year=pred.weights$year, source='Prediction'))
w.comp=w.comp[w.comp$w<=quantile(w.comp$w, probs = 0.995),]

p.comb=ggplot(data=w.comp)+
  geom_density(aes(x=(w), group=source, fill=source), alpha=0.2)+
  scale_fill_viridis_d()+
  ylab('Density')+
  xlab('Weight (grams)')+
  theme(legend.position = 'bottom')+
  labs(fill='Data');p.comb

mean.preds=pred.weights%>%
  dplyr::group_by(year)%>%
  dplyr::summarise(w.mu=mean(log(M7)), w.sd=sd(log(M7)), source='Pred')
mean.obs=obs.long%>%
  dplyr::group_by(year)%>%
  dplyr::summarise(w.mu=mean(log(weig_g)), w.sd=sd(log(weig_g)), source='Obs')

p4=rbind(mean.preds, mean.obs)%>%
  ggplot(aes(x=year, y=w.mu, color=source, group=source))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=w.mu-w.sd,ymax=w.mu+w.sd, fill=source), color=NA,alpha=0.1)+
  theme(legend.position = 'none')+
  scale_color_viridis_d()+
  ylab('log Weight (g)')+
  xlab('Year')+
  scale_fill_viridis_d();p4

p.temp=obs.temp%>%
  dplyr::mutate(season=ifelse(month%in%5:6, 'Spring',ifelse(month%in%7:9,'Summer','Autumn')))%>%
  dplyr::group_by(season, year)%>%
  dplyr::summarise(mu.t=mean(mu.t))%>%
  ggplot(aes(x=year, y=mu.t, color=season))+
  geom_line()+
  geom_point()+
  scale_color_viridis_d()+
  theme(legend.position = 'bottom')+
  labs(fill='Season', color='Season')+
  ylab('Mean seasonal T°')+
  xlab('Year');p.temp

p.vars=ggpubr::ggarrange(p4,p.temp, ncol=1, labels = c("b)", "c)"));p.vars
p.fin=ggpubr::ggarrange(p.comb, p.vars, labels = c("a)"), common.legend = T)
ggsave('results/plots/Figure6.jpeg', width=20, height=15, units='cm', dpi=500)


# save results to reproduce plots

dat7d=obs.temp%>%
  dplyr::mutate(season=ifelse(month%in%5:6, 'Spring',ifelse(month%in%7:9,'Summer','Autumn')))%>%
  dplyr::group_by(season, year)%>%
  dplyr::summarise(mu.t=mean(mu.t)) # ptem
dat7b=w.comp # pobs
dat7c=rbind(mean.preds, mean.obs) # p4
saveRDS(list(d7b=dat7b, d7c=dat7c, d7d=dat7d), file='results/plots/plotdata.RDS')


