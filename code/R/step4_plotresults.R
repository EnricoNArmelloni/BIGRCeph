# Description ####
# this script extract posterior MCMC draws, provide summary results and generate the images accompaining the manuscript
rm(list=ls())

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
read.dir='results'
model.dir='results/model_fits'
save.dir='results'
setwd(project.dir)

# libraries and functions
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
source('code/R/supporting_functions.R')

### get data ####
best.M1=readRDS(file.path(model.dir, 'M1_4_J_fit.RDS'))
best.M2=readRDS(file.path(model.dir, 'M2_4_Z_fit.RDS'))
M.data=read_csv("data/input_data.csv")
M1.data=read_csv("data/M1_model_data.csv")
M2.data=read_csv("data/M2_model_data.csv")
M12.data=rbind(M1.data, M2.data)

### get posteriors ####
M1.posteriors=format.posteriors(best.M1)
M2.posteriors=format.posteriors(best.M2)

M1.post.tab=data.frame(mean=apply(M1.posteriors,2,mean))
M1.post.tab$sd=apply(M1.posteriors,2,sd)
M1.post.tab$var=rownames(M1.post.tab)

M1.post.tab$distr=paste0(round(M1.post.tab$mean, digits=2),' , ', round(M1.post.tab$sd, digits=2) )
write.csv(M1.post.tab, 'results/other_tables/m1_post_avg.csv', row.names = F)

M2.post.tab=data.frame(mean=apply(M2.posteriors,2,mean))
M2.post.tab$sd=apply(M2.posteriors,2,sd)
M2.post.tab$var=rownames(M2.post.tab)
M2.post.tab$distr=paste0('~ N(', round(M2.post.tab$mean, digits=2),',', round(M2.post.tab$sd, digits=2) ,')')
write.csv(M2.post.tab, 'results/other_tables/m2_post_avg.csv', row.names = F)

# Conditional effects ####
R.vector=seq(0,1,0.05)
T.vector=data.frame(t.transf=seq(min(M12.data$z2), max(M12.data$z2), length.out=1000),
                       temp=seq(min(M12.data$Temperature), max(M12.data$Temperature), length.out=1000))%>%
  dplyr::mutate(xvar2=round(t.transf, digits=2))%>%
  distinct(xvar2, .keep_all = T)


## Ration
R.effect.M1=NULL
for(i in 1:length(R.vector)){
  i.res=data.frame(xvar=R.vector[i], 
                   resp=M1(a=M1.posteriors$a,
                           b1=M1.posteriors$bF,
                           b2=M1.posteriors$bT1,
                           b3=M1.posteriors$bT2,
                           sigma=M1.posteriors$`sigma[28]`,
                           mu=M1.posteriors$mu,
                           t=16,
                           tem=0.425,
                           Fi=M1.posteriors$`F[1]`,
                           R=R.vector[i]))
  R.effect.M1=rbind(R.effect.M1, i.res)
}
R.effect.M1=R.effect.M1%>%
  dplyr::group_by(xvar)%>%
  dplyr::summarise(resp.mu=quantile(resp, probs=c(0.5)),
                   resp.lo=quantile(resp, probs=c(0.05)),
                   resp.hi=quantile(resp, probs=c(0.95)))

p.r.m1=base.plot.curve(R.effect.M1)+
  xlab("Ration / Ad libitum")+
  ylab(expression(paste(IGR['t = 16'], ' (% bw *',  day^-1,")")));p.r.m1

R.effect.M2=NULL
for(i in 1:length(R.vector)){
  i.res=data.frame(xvar=R.vector[i], 
                   resp=M2(A=M2.posteriors$A,
                           b1=M2.posteriors$bF,
                           b4a=M2.posteriors$bTA,
                           b2=M2.posteriors$bTB1,
                           b3=M2.posteriors$bTB2,
                           B=M2.posteriors$`B[39]`,
                           C=M2.posteriors$C,
                           t=100,
                           tem=0.425,
                           Fi=M2.posteriors$`F[1]`,
                           R=R.vector[i]))
  R.effect.M2=rbind(R.effect.M2, i.res)
}

R.effect.M2=R.effect.M2%>%
  dplyr::group_by(xvar)%>%
  dplyr::summarise(resp.mu=quantile(resp, probs=c(0.5)),
                   resp.lo=quantile(resp, probs=c(0.05)),
                   resp.hi=quantile(resp, probs=c(0.95)))

p.r.m2=base.plot.curve(R.effect.M2)+
  xlab("Ration / Ad libitum")+
  ylab(expression(paste(IGR['t = 50'], ' (% bw *',  day^-1,")")))

## Temperature
T.effect.M1=NULL
for(i in 1:nrow(T.vector)){
  i.res=data.frame(xvar=T.vector[i,]$temp, 
                   resp=M1(a=M1.posteriors$a,
                           b1=M1.posteriors$bF,
                           b2=M1.posteriors$bT1,
                           b3=M1.posteriors$bT2,
                           sigma=M1.posteriors$`sigma[28]`,
                           mu=M1.posteriors$mu,
                           t=16,
                           tem=T.vector[i,]$xvar2,
                           Fi=M1.posteriors$`F[1]`,
                           R=0.5))
  T.effect.M1=rbind(T.effect.M1, i.res)
}
T.effect.M1=T.effect.M1%>%
  dplyr::group_by(xvar)%>%
  dplyr::summarise(resp.mu=quantile(resp, probs=c(0.5)),
                   resp.lo=quantile(resp, probs=c(0.05)),
                   resp.hi=quantile(resp, probs=c(0.95)))
p.t.m1=base.plot.curve(T.effect.M1)+
  xlab("T°")+
  ylab(expression(paste(IGR['t = 16'], ' (% bw *',  day^-1,")")))

T.effect.M2=NULL
for(i in 1:nrow(T.vector)){
  i.res=data.frame(xvar=T.vector[i,]$temp, 
                   A=M2.posteriors$A*(T.vector[i,]$xvar2+0.5)^M2.posteriors$bTA*(0.5*M2.posteriors$`F[1]`)^M2.posteriors$bF,
                   B=M2.posteriors$`B[39]`*((T.vector[i,]$xvar2+0.5)*M2.posteriors$bTB1 + M2.posteriors$bTB2*(T.vector[i,]$xvar2+0.5)^2 ))
  
  T.effect.M2=rbind(T.effect.M2, i.res)
}

T.effect.M2=T.effect.M2%>%
  dplyr::mutate(A=A/max(A),B=B/max(B))%>%
  pivot_longer(-xvar)%>%
  dplyr::group_by(xvar,name)%>%
  dplyr::summarise(resp.mu=quantile(value, probs=c(0.5)),
                   resp.lo=quantile(value, probs=c(0.05)),
                   resp.hi=quantile(value, probs=c(0.95)))

p.t.m2=ggplot(data=T.effect.M2, aes(x=xvar, y=resp.mu, group=name))+
  geom_line()+
  geom_ribbon(aes(ymin=resp.lo, ymax=resp.hi, fill=name), alpha=0.2)+
  theme_minimal()+
  xlab("T°")+
  #ylab(expression(paste(IGR['t = 50'], ' (% bw *',  day^-1,")")))+
  ylab('Relative effect of parameter')+
  scale_color_viridis_d(drop = FALSE)+
  scale_fill_viridis_d(drop = FALSE)+
  theme(legend.position = 'bottom')+
  annotate(x = 25, y = 0.4, geom = "text", label = expression(italic(b['TB'])))+
  annotate(x = 25, y = 0.1, geom = "text", label = expression(italic(b['TA'])));p.t.m2

# Prey type
base.plot.food=function(xdat){
  ggplot(data=xdat, aes(x=post, fill=diet))+
    geom_density( alpha=0.7)+
    theme_minimal()+
    theme(legend.position = 'bottom',
          legend.direction = "horizontal")+
    scale_fill_viridis_d(drop = FALSE)+
    guides(fill = guide_legend(nrow = 2))
  
}

index.F.M1=distinct(M1.data, Index.food, main_prey)
M1.post.F=M1.posteriors[,grep('F', names(M1.posteriors))]#
M1.post.F=M1.post.F[,-grep('bF', names(M1.post.F))]#
names(M1.post.F)=index.F.M1$main_prey
M1.post.F=M1.post.F%>%
  dplyr::mutate(id=seq(1:nrow(.)))%>%
  pivot_longer(-id, names_to = 'diet', values_to = 'post')

index.F.M2=distinct(M2.data, Index.food, main_prey)
M2.post.F=M2.posteriors[,grep('F', names(M2.posteriors))]#
M2.post.F=M2.post.F[,-grep('bF', names(M2.post.F))]#
names(M2.post.F)=index.F.M2$main_prey
M2.post.F=M2.post.F%>%
  dplyr::mutate(id=seq(1:nrow(.)))%>%
  pivot_longer(-id, names_to = 'diet', values_to = 'post')

M2.post.F%>%dplyr::group_by(diet)%>%dplyr::summarise(mean(post))

uniq.diets=unique(c(unique(M1.post.F$diet), unique(M2.post.F$diet)))
my.colors=viridis::viridis(n=length(uniq.diets))
my.color.2=data.frame(prey=uniq.diets,
                      color=my.colors)

p.f.m1=base.plot.food(xdat=M1.post.F)+
  scale_fill_manual(  values = my.color.2[my.color.2$prey %in% M1.post.F$diet,]$color,
                      breaks = my.color.2[my.color.2$prey %in% M1.post.F$diet,]$prey)+
  xlab('Posterior value')

p.f.m2=base.plot.food(xdat=M2.post.F)+
  scale_fill_manual(  values = my.color.2[my.color.2$prey %in% M2.post.F$diet,]$color,
                      breaks = my.color.2[my.color.2$prey %in% M2.post.F$diet,]$prey)+
  xlab('Posterior value')+
  labs(fill='Prey Type')

# combine plots
plot.M1=ggpubr::ggarrange(p.r.m1,p.t.m1,p.f.m1, ncol=3, legend = F, labels = c("a)", "b)", "c)"));plot.M1
plot.M2=ggpubr::ggarrange(p.r.m2,p.t.m2,p.f.m2, ncol=3, legend = F, labels = c("d)", "e)", "f)"));plot.M2
plot.effect=ggpubr::ggarrange(plot.M1, plot.M2,ncol=1, 
                              common.legend = T, legend.grob = ggpubr::get_legend(p.f.m2),
                      legend='bottom');plot.effect
ggsave(plot=plot.effect, 'results/plots/covars_effect.jpeg', width = 20, height = 15, units='cm')

## IGR curve ####
M12.data$temp=M12.data$Temperature
M12.data$ration=M12.data$Rcomma
unique.temp=M12.data%>%distinct(temp,z2)%>%
  dplyr::filter(temp%in%c(12.5,19,27))
unique.ratio=M12.data%>%distinct(ration,z1)%>%
  arrange(ration)
unique.ratio=unique.ratio[c(1,17,23),]
pred.grid=expand.grid(temp=unique(unique.temp$temp), ration=unique(unique.ratio$ration))
pred.grid=left_join(pred.grid, unique.temp, by=c('temp'))%>%
  left_join(.,unique.ratio, by='ration')

x.dat=M.data[M.data$main_prey=='medium_shrimp',]
x.dat$IGR.mu=x.dat$igr.mu
x.dat$R=x.dat$ration
store.pred.M1=store.pred.M2=store.pred.comb=list()
time.vec=seq(0,120,5)
for(j in 1:nrow(pred.grid)){
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
                              tem=pred.grid[j,]$z2,
                              Fi=M1.posteriors$`F[4]`,
                              R=pred.grid[j,]$z1))
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
                               tem=pred.grid[j,]$z2,
                               Fi=M2.posteriors$`F[2]`,
                               R=pred.grid[j,]$z1))
    m2.pred$pred=as.numeric(mapply(function(mu, sigma) rnorm(n = 1, mean = mu, sd = sigma), m2.pred$pred, 
                                   apply(M2.posteriors[,grep('obs_erro', names(M2.posteriors))],1,mean)))
    m2.store[,i]=m2.pred$pred
    
  }
  m1.store=as.data.frame(m1.store)
  names(m1.store)=paste0('t', time.vec)
  m2.store=as.data.frame(m2.store)
  names(m2.store)=paste0('t', time.vec)
  m.comb=as.data.frame(apply_sigmoid_row(m1.store,m2.store,time.vec))
  m.comb$z1=m1.store$z1=m2.store$z1=pred.grid[j,]$z1
  m.comb$z2=m1.store$z2=m2.store$z2=pred.grid[j,]$z2
  m.comb$R=m1.store$R=m2.store$R=pred.grid[j,]$ration
  m.comb$temp=m1.store$temp=m2.store$temp=pred.grid[j,]$temp
  store.pred.comb[[j]]=m.comb
  store.pred.M1[[j]]=m1.store
  store.pred.M2[[j]]=m2.store
}

M1.resp=plyr::ldply(store.pred.M1)%>%
  pivot_longer(-c(R, temp, z1,z2), names_to = 't', values_to = 'IGR')%>%
  dplyr::mutate(t=as.numeric(str_remove(t,'t')))%>%
  dplyr::filter(t<44)%>%
  dplyr::group_by(t,temp,R,z1,z2)%>%
  dplyr::summarise(IGR.mu=mean(IGR, probs=c(0.5)),
                   IGR.lo=quantile(IGR, probs=c(0.1)),
                   IGR.hi=quantile(IGR, probs=c(0.9)))

M2.resp=plyr::ldply(store.pred.M2)%>%
  pivot_longer(-c(R, temp, z1,z2), names_to = 't', values_to = 'IGR')%>%
  dplyr::mutate(t=as.numeric(str_remove(t,'t')))%>%
  dplyr::filter(t>16)%>%
  dplyr::group_by(t,temp,R,z1,z2)%>%
  dplyr::summarise(IGR.mu=mean(IGR, probs=c(0.5)),
                   IGR.lo=quantile(IGR, probs=c(0.1)),
                   IGR.hi=quantile(IGR, probs=c(0.9)))

M12.resp=rbind(M1.resp%>%dplyr::mutate(model='M1'),
          M2.resp%>%dplyr::mutate(model='M2'))

Mcomb.resp=plyr::ldply(store.pred.comb)%>%
  pivot_longer(-c(R, temp, z1,z2), names_to = 't', values_to = 'IGR')%>%
  dplyr::mutate(t=as.numeric(str_remove(t,'t')))%>%
  dplyr::group_by(t,temp,R,z1,z2)%>%
  dplyr::summarise(IGR.mu=mean(IGR, probs=c(0.5)),
                   IGR.lo=quantile(IGR, probs=c(0.1)),
                   IGR.hi=quantile(IGR, probs=c(0.9)))

baseplot.response=function(xdat){
 if('model'%in% names(xdat)){
   ggplot(data=xdat, aes(x=t, y=IGR.mu, group=model))+
     geom_line()+
     facet_grid(rows = vars(temp), cols=vars(round(R, digits=2)))+
     
     geom_ribbon(aes(ymin=IGR.lo, ymax=IGR.hi, fill=model),
                 alpha=0.2)+
     theme_minimal()+
     ylab('IGR (%bw*day-1')+
     theme(legend.position = 'bottom')+
     scale_color_viridis_d()+
     scale_fill_viridis_d()
 }else{
   ggplot(data=xdat, aes(x=t, y=IGR.mu))+
     geom_line()+
     facet_grid(rows = vars(temp), cols=vars(round(R, digits=2)))+
     
     geom_ribbon(aes(ymin=IGR.lo, ymax=IGR.hi),
                 alpha=0.2)+
     theme_minimal()+
     ylab('IGR (%bw*day-1')+
     theme(legend.position = 'bottom')+
     scale_color_viridis_d()+
     scale_fill_viridis_d() 
 } 
}

pc1=baseplot.response(M12.resp)+
  #geom_vline(xintercept = 16, linetype=2)+
  labs(fill='Model', color='Model')+
  xlab('Time (DAH)')+
  scale_x_continuous(breaks = seq(0,120,30))+
  ylab(expression(paste(IGR, ' (% body weight *',  day^-1,")")));pc1
tgr=textGrob("R'", gp = gpar(fontsize = 8, fontface = "bold"))
tgc=textGrob("T°", gp = gpar(fontsize = 8, fontface = "bold"), rot=-90)
pc1=grid.arrange(tgr,
                 arrangeGrob(pc1,tgc,  nrow = 1, widths  = c(1, 0.05)),  # Column label at the top
                 
                 ncol = 1, heights = c(0.05, 1))  # Row label on the left

pc2=baseplot.response(Mcomb.resp)+
  labs(fill='Model', color='Model')+
  xlab('Time (DAH)')+
  scale_x_continuous(breaks = seq(0,120,30))+
  #geom_point(data=x.dat, alpha=0.1, color='black')+
  ylab(expression(paste(IGR, ' (% body weight *',  day^-1,")")));pc2

pc3=ggpubr::ggarrange(pc1,pc2, ncol=2, labels = c("a)", "b)"))
#ggsave(plot=pc3, 'results/plots/prediction_combined_pubr.jpeg', width = 20, height = 10, units='cm')


p4=ggplot(data=Mcomb.resp[round(Mcomb.resp$R, digits=2)==0.86,], aes(x=t, y=IGR.mu))+
  geom_line(aes(color=factor(temp)))+
  geom_ribbon(aes(ymin=IGR.lo, ymax=IGR.hi,fill=factor(temp)),
              alpha=0.2)+
  theme_minimal()+
  ylab('IGR (%bw*day-1')+
  theme(legend.position = 'bottom')+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab('Time (DAH)')+
  labs(color='T (C°)', fill='T (C°)')+
  scale_x_continuous(breaks = seq(0,120,30))
pc4=ggpubr::ggarrange(pc1,p4, ncol=2, labels = c("a)", "b)"))
ggsave(plot=pc4, 'results/plots/prediction_combined_pubr2.jpeg', width = 20, height = 10, units='cm')


### validation beyond
alldat= read_csv("data/input_data.csv")
x.dat=alldat[alldat$t>=120,]
x.dat[is.na(x.dat$Rcomma),]$Rcomma=1

#x.dat=x.dat[x.dat$main_prey!='crayfish',]
extra.exp=unique(x.dat$exp_id)
x.dat$IGR.mu=x.dat$igr.mu

conv2=M12.data%>%distinct(temp,z2)%>%
  dplyr::filter(temp %in% unique(x.dat$Temperature))
conv1=M12.data%>%distinct(ration,z1)%>%
  dplyr::filter(ration %in% unique(x.dat$Rcomma))
x.dat=x.dat[x.dat$exp_id!='syk_06_3_3',]
store.pred.M2=list()

extra.exp=extra.exp[extra.exp!='syk_06_3_3']
for(j in 1:length(extra.exp)){
  j.exp=extra.exp[j]
  j.dat=x.dat[x.dat$exp_id==j.exp,]
  
  j.temp=conv2[conv2$temp==unique(j.dat$Temperature),]$z2
  j.R=conv1[conv1$ration==unique(j.dat$Rcomma),]$z1
  j.diet=unique(j.dat$main_prey)
  j.diet=M2.post.F[M2.post.F$diet==j.diet,]$post
  time.vec=unique(j.dat$t)
 
  m2.store=matrix(nrow = nrow(M2.posteriors), ncol=length(time.vec))
  for(i in 1:length(time.vec)){  
    m2.pred=data.frame(pred=M2(A=M2.posteriors$A,
                               b1=M2.posteriors$bF,
                               b4a=M2.posteriors$bTA,
                               b2=M2.posteriors$bTB1,
                               b3=M2.posteriors$bTB2,
                               B=M2.posteriors$`B[39]`,
                               C=M2.posteriors$C,
                               t=time.vec[i],
                               tem=j.temp,
                               Fi=j.diet,
                               R=j.R))
    m2.pred$pred=as.numeric(mapply(function(mu, sigma) rnorm(n = 1, mean = mu, sd = sigma), m2.pred$pred,
                                   apply(M2.posteriors[,grep('obs_erro', names(M2.posteriors))],1,mean)))
    m2.store[,i]=m2.pred$pred
    
  }
  m2.store=as.data.frame(m2.store)
  names(m2.store)=paste0('t', time.vec)
  m2.store$exp_id=j.exp
  m2.store$Temperature=unique(j.dat$Temperature)
  m2.store$Rcomma=unique(j.dat$Rcomma)
  store.pred.M2[[j]]=m2.store
}

M2.resp=plyr::ldply(store.pred.M2)%>%
  pivot_longer(-c(exp_id, Temperature, Rcomma), names_to = 't', values_to = 'IGR')%>%
  dplyr::filter(!is.na(IGR))%>%
  dplyr::mutate(t=as.numeric(str_remove(t,'t')))%>%
  dplyr::group_by(t,Temperature, Rcomma)%>%
  dplyr::summarise(IGR.mu=mean(IGR, probs=c(0.5)),
                   IGR.lo=quantile(IGR, probs=c(0.1)),
                   IGR.hi=quantile(IGR, probs=c(0.9)))

xx.dat=left_join(x.dat, M2.resp, by=c('Temperature', 'Rcomma', 't'))
xx.dat$residuals=xx.dat$igr.mu-xx.dat$IGR.mu.y



pcc=ggplot(data=xx.dat, aes(group=1))+
  #facet_wrap(~exp_id)+
  facet_grid(rows = vars(Temperature), cols=vars(round(Rcomma, digits=2)))+
  geom_point(aes(x=t, y=igr.mu))+
  geom_line(aes(x=t, y=IGR.mu.y, group = 1))+
  theme(legend.position = 'bottom')+
  scale_color_viridis_d()+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  xlab('Time')+
  labs(fill= 'Intake (relative)')+
  xlab('Time (DAH)')+
  labs(color='Temperature (C°)')+
  geom_ribbon(data=M2.resp,aes(x=t, ymin=IGR.lo, ymax=IGR.hi, group = 1), alpha=0.2)+
  ylab(expression(paste(IGR, ' (% body weight *',  day^-1,")")));pcc


tgr=textGrob("R'", gp = gpar(fontsize = 8, fontface = "bold"))
tgc=textGrob("T°", gp = gpar(fontsize = 8, fontface = "bold"), rot=-90)
pcc=grid.arrange(tgr,
                 arrangeGrob(pcc,tgc,  nrow = 1, widths  = c(1, 0.05)),  # Column label at the top
                 
                 ncol = 1, heights = c(0.05, 1));pcc  # Row label on the left

ggsave(plot=pcc, 'results/plots/prediction_beyond.jpeg', width = 10, height = 10, units='cm')

ej=mean(abs(xx.dat$residuals))
MAD=mean(abs(xx.dat$igr.mu-mean(xx.dat$igr.mu)))
ej/MAD ### MASE

# Compute Coverage Probability
df <- xx.dat %>%
  mutate(inside_CI = ifelse(igr.mu >= IGR.lo & igr.mu <= IGR.hi, 1, 0)) 

coverage_prob <- mean(df$inside_CI, na.rm = TRUE)
cat("Coverage Probability:", round(coverage_prob, 3), "\n")

# other plots ####

resp.type1=function(alpha,beta,x){
  alpha*x^beta
}
x=seq(0,1, 0.005)
resp1.df=expand.grid(x=x, beta=c(0.5,1,1.5),alpha=c(1))%>%
  dplyr::mutate(y=resp.type1(alpha=alpha,beta=beta,x=x))
p1=ggplot(data=resp1.df, aes(x=x,y=y, linetype=factor(beta)))+
  geom_line()+
  theme_minimal()+
  ylab('Covariate effect')+
  xlab('Covariate value')+
  labs(linetype=expression(italic('b')))+
  theme(legend.position = 'bottom')

#ggsave('results/plots/effect1.jpeg', width =8, height = 8, units='cm')
resp.type2=function(alpha,beta1,beta2,x){
  alpha+beta1*x+beta2*x^2
}
x=seq(0,1, 0.005)
resp2.df=expand.grid(x=x, beta1=c(0.5,1,1.5),beta2=c(0,-0.5,-1),alpha=c(1))%>%
  dplyr::mutate(y=resp.type2(alpha=alpha,beta1=beta1,beta2=beta2,x=x))
p2=ggplot(data=resp2.df, aes(x=x,y=y, color=factor(beta2), linetype=factor(beta1)))+
  geom_line()+
  #facet_wrap(~beta1)+
  theme_minimal()+
  ylab('Covariate effect')+
  xlab('Covariate value')+
  theme(legend.position = 'bottom')+
  labs(linetype=expression(italic('b'[1])),
       color=expression(italic('b'[2])))+
  scale_color_viridis_d()

p3=ggpubr::ggarrange(p1,p2, labels = c("a)", "b)"))

ggsave(plot=p3,'results/plots/effect2.jpeg', width = 25, height = 8, units='cm')






