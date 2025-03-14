rm(list=ls())
# Description ####
# this script applies model diagnostics to stan fits. Most of the functions applied are custom functions (largely based on bayesplot, loo and ggplot), which code is stored in the supporting_functions.R file

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
read.dir='results/model_fits'
save.dir='results/diags'
setwd(project.dir)

# libraries and functions
library(tidyverse)
library(loo)
library(bayesplot)
library(readxl)
source('code/R/supporting_functions.R')

# options
color_scheme_set("viridis")
theme_set(theme_bw())

# M1 ####
M1.dat=read.csv(file.path('data', 'M1_model_data.csv'))
model.list=list.files(path = read.dir, pattern = '.RDS')
model.list=model.list[grep('M1',model.list)]

prior.M1 <- data.frame(
  pr_mean_mu=(20), pr_mean_mu_sd=(1.5),
  pr_mean_sigma= (15), pr_mean_sigma_sd=(10),
  # prior for mu and sigma when hierarchical
  pr_hyper_mu=20, pr_hyper_mu_sd=2.5,
  pr_hyper_mu_deviation=1.5,pr_hyper_mu_deviation_sd=0.5,
  pr_hyper_sigma=15,   pr_hyper_sigma_sd=10,
  pr_hyper_sigma_deviation=1.5,pr_hyper_sigma_deviation_sd=0.5,
  # priors for covariates
  pr_a=(1), pr_a_sd= exp(0.5), # updated
  pr_bpow=(1), pr_bpow_sd=exp(0.5),
  pr_bquad=(1),pr_bquad_sd=exp(0.5),
  pr_diet=0, pr_diet_sd=1, 
  # prior for observation error
  pr_sigma_obs=5, pr_sigma_obs_sd=3)

### launch diags
for(x.models in 1:length(model.list)){
  
  # load model
  fit_mcmc=readRDS(file.path(read.dir, model.list[x.models]))
  
  # set up for labels and directory
  model.lab=str_remove(model.list[x.models],'_fit.RDS')
  loop.savedir=file.path(save.dir, model.lab)
  dir.create(loop.savedir)
  loop.savedir=paste0(loop.savedir,'/')
  hierarchy=substr(model.lab, nchar(model.lab),nchar(model.lab))
  
  # get quantities from the model
  x.summary=fit_mcmc$summary()
  posterior_draws = fit_mcmc$draws()
  lp_cp=log_posterior(fit_mcmc)
  np_cp <- nuts_params(fit_mcmc)
  
  # check the parameters
  x.pars=unique(gsub("\\[.*?\\]", "", x.summary$variable))
  if(hierarchy=='J'){
    fixed.pars=x.pars[which(x.pars %in% c('a','bF','bT','bT1','bT2','sigma_obs','mu'))]
    random.pars=c("mean_sigma","sigma_sigma", paste0('sigma[',1:28,']'))
  }else if(hierarchy=='I'){
    fixed.pars=x.pars[which(x.pars %in% c('a','bF','bT','bT1','bT2','sigma_obs','sigma'))]
    random.pars=c("mean_mu","sigma_mu", paste0('mu[',1:28,']'))
  }else if(hierarchy=='N'){
    fixed.pars=x.pars[which(x.pars %in% c('a','bF','bT','bT1','bT2','sigma_obs','sigma','mu'))]
    random.pars=NA
  }
  
  # mcmc diags
  mc.diags=mcmc.diags(x.fit = fit_mcmc,
                 label=model.lab,
                 x.post=posterior_draws,
                 np_cp = np_cp, lp_cp=lp_cp, fixed.pars = fixed.pars, random.pars = random.pars,
                 n.iter = nrow(lp_cp))
  
  # priors vs posteriors
  label_list=list(expression(sigma[obs]),
                  expression(mu),
                  'a',
                  expression(b[F]),
                  expression(b[T1]),
                  expression(b[T2]),  expression(Hyper~sigma[mu]), expression(Hyper~sigma[sigma]))
  
  posterior.diags(fit_mcmc=fit_mcmc,
                      prior.list = prior.M1,
                      label=model.lab,
                      fixed.pars = fixed.pars,
                      random.pars = random.pars,
                  label_list = label_list,
                      save.dir=save.dir)
  
 fi.diags=rmse.fit(posterior_draws = posterior_draws,
           y.dat=M1.dat,
           x.fit=fit_mcmc,
           label=model.lab, 
           save.dir=loop.savedir)
 
 ## build results row
 x.res=data.frame(estimate=t(mc.diags), se=NA, par=names(mc.diags))
 diags.df=rbind(x.res, fi.diags)
 diags.df$model=model.lab
 write.csv(diags.df, paste0(loop.savedir, 'restable.csv'), row.names = F)
}

# performance metrics
M1.models=list.files(save.dir, pattern='M1')
store=NULL
for(i in 1:length(M1.models)){
  x.res=read.csv(file.path(save.dir, M1.models[i], 'restable.csv'))
  x.res$model=M1.models[i]
  store=rbind(store, x.res)
}
store.tab=store%>%
  dplyr::select(-se)%>%
  pivot_wider(names_from = par, values_from = estimate)
#write.csv(store.tab, 'results/metricsM1.csv', row.names = F)

# get posterior of the best model
best.model='M1_4_J'
best.fit=readRDS(file.path(read.dir, model.list[grep(best.model, model.list)]))
best.posteriors = best.fit$draws()
best.mu=best.posteriors[,,'mu']
#best.mu=best.posteriors[,,'mu[28]']
mu.estimation=ceiling(mean(unlist(data.frame(best.mu))))


# M2 ####
M2.dat=read.csv(file.path('data', 'M2_model_data.csv'))
model.list=list.files(path = read.dir, pattern = '.RDS')
model.list=model.list[grep('M2',model.list)]

prior.M2 <- data.frame(
                      # prior for mu and sigma when is shared
                      pr_mean_A=10, pr_mean_A_sd=10,
                      pr_mean_B= 10, pr_mean_B_sd=10,
                      # prior for mu and sigma when hierarchical
                      pr_hyper_A=5, pr_hyper_A_sd=1.5,
                      pr_hyper_A_deviation=1.5,pr_hyper_A_deviation_sd=1,
                      pr_hyper_B=10,   pr_hyper_B_sd=1.5,
                      pr_hyper_B_deviation=1.5,pr_hyper_B_deviation_sd=1,
                      # priors for covariates
                      pr_C=(0.05), pr_C_sd= exp(0.5),
                      pr_bpow=(1), pr_bpow_sd=(2),
                      pr_bquad=(3),pr_bquad_sd=(2),
                      pr_diet=log(0.5), pr_diet_sd=0.5, 
                      # prior for observation error
                      pr_sigma_obs=5, pr_sigma_obs_sd=3)

### launch diags
for(x.models in 1:length(model.list)){
  
  # load model
  fit_mcmc=readRDS(file.path(read.dir, model.list[x.models]))
  # set up for labels and directory
  model.lab=str_remove(model.list[x.models],'_fit.RDS')
  loop.savedir=file.path(save.dir, model.lab)
  dir.create(loop.savedir)
  loop.savedir=paste0(loop.savedir,'/')
  hierarchy=substr(model.lab, nchar(model.lab),nchar(model.lab))
  
  # get quantities from the model
  x.summary=fit_mcmc$summary()
  posterior_draws = fit_mcmc$draws()
  lp_cp=log_posterior(fit_mcmc)
  np_cp <- nuts_params(fit_mcmc)
  
  # check the parameters
  x.pars=unique(gsub("\\[.*?\\]", "", x.summary$variable))
  
  if(hierarchy=='Z'){
    fixed.pars=x.pars[which(x.pars %in% c('A', 'C', 'bF','bTA','bTB','bTB1','bTB2', 'sigma_obs','mu'))]
    random.pars=x.pars[which(x.pars %in% c("mean_B","sigma_B",'B'))]
    random.pars=c("mean_B","sigma_B", paste0('B[',1:39,']'))
  }else if(hierarchy=='K'){
    fixed.pars=x.pars[which(x.pars %in% c('B', 'C', 'bF','bTA','bTB','bTB1','bTB2','sigma_obs','mu'))]
    random.pars=x.pars[which(x.pars %in% c("mean_A","sigma_A",'A'))]
  }else if(hierarchy=='N'){
    fixed.pars=x.pars[which(x.pars %in%c('A', 'B', 'C', 'bF','bTA','bTB','bTB1','bTB2','sigma_obs','mu'))]
    random.pars=NA
  }
  
  # mcmc diags
  mc.diags=mcmc.diags(x.fit = fit_mcmc,
                          label=model.lab,
                          x.post=posterior_draws,
                          np_cp = np_cp, lp_cp=lp_cp, fixed.pars = fixed.pars, random.pars = random.pars,
                          n.iter = nrow(lp_cp))
  
  # priors vs posteriors
  label_list=list('A',
                  expression(sigma[obs]),
                  expression(b[TA]),
                  'C',
                  expression(b[F]),
                  expression(b[TB1]),
                  expression(b[TB2]),  expression(Hyper~B[mu]), expression(Hyper~B[sigma]))
  
  posterior.diags(fit_mcmc = fit_mcmc,
                      prior.list = prior.M2,
                      label=model.lab,
                      fixed.pars = fixed.pars,
                  label_list = label_list,
                      random.pars = random.pars,
                      save.dir=save.dir)
  
  fi.diags=rmse.fit(posterior_draws = posterior_draws,
                    y.dat=M2.dat,
                    x.fit=fit_mcmc,
                    label=model.lab, 
                    save.dir=loop.savedir)
  
  ## build results row
  x.res=data.frame(estimate=t(mc.diags), se=NA, par=names(mc.diags))
  diags.df=rbind(x.res, fi.diags)
  diags.df$model=model.lab
  write.csv(diags.df, paste0(loop.savedir, 'restable.csv'), row.names = F)
}

# performance metrics
M2.models=list.files(save.dir, pattern='M2')
store=NULL
for(i in 1:length(M2.models)){
  x.res=read.csv(file.path(save.dir, M2.models[i], 'restable.csv'))
  x.res$model=M2.models[i]
  store=rbind(store, x.res)
}
store.tab=store%>%
  dplyr::select(-se)%>%
  pivot_wider(names_from = par, values_from = estimate)
write.csv(store.tab, 'results/metricsM2.csv', row.names = F)

