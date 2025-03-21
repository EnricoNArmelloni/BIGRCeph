# Description ####
# Apply model variants (.stan models) and save all the outputs in a separate folder (outside the project, for space constrains). In this version only the best model variants are activated for space reasons. To apply other models, just remove the # simbol in lines 90-103 and 205-219 
rm(list=ls())

# set paths
project.dir="C:/Users/e.armelloni/OneDrive/Lezioni/Lavoro/PhD/Activities/Task5/analysis/BIGRCeph"
data.dir='data'
models.dir='code/stan'
save.dir='results/model_fits'
setwd(project.dir)

# libraries and functions
library(cmdstanr)

# options
options(mc.cores = parallel::detectCores())

# read data
x.dat= read_csv(file.path(data.dir, "input_data.csv"))
x.dat=x.dat[x.dat$t<=120,]
x.dat=x.dat[!is.na(x.dat$Rcomma),]

# min-max transformation 
x.dat$z1=((x.dat$Rcomma)-min(x.dat$Rcomma))/(max(x.dat$Rcomma)-min(x.dat$Rcomma))
x.dat$z2=((x.dat$Temperature)-min(x.dat$Temperature))/(max(x.dat$Temperature)-min(x.dat$Temperature))


# M1 models ####
y.dat.j=x.dat[x.dat$t<40,]
set.seed(14)
y.dat.j=y.dat.j[sample(nrow(y.dat.j), nrow(y.dat.j)),]

# create indices for hierarchy ####
y.dat=y.dat.j
x.id.t=unique(y.dat$exp_id)
y.dat$Index.t=NA
for(i in 1:nrow(y.dat)){
  i.id=y.dat[i,]$exp_id
  y.dat[i,]$Index.t=which(x.id.t==i.id)
}

x.id.e=unique(paste0(y.dat$experiment))
y.dat$Index.e=NA
for(i in 1:nrow(y.dat)){
  i.id=paste0(y.dat[i,]$experiment)
  y.dat[i,]$Index.e=which(x.id.e==i.id)
}
y.dat$Index.no=1

# create indices for food items ####
x.id.food=unique(y.dat$main_prey)
y.dat$Index.food=NA
for(i in 1:nrow(y.dat)){
  i.id=y.dat[i,]$main_prey
  y.dat[i,]$Index.food=which(x.id.food==i.id)
}

write.csv(y.dat, file.path(data.dir, 'M1_model_data.csv'), row.names = F)

# stan data
# priors used in the paper ####
data.M1= list(N = nrow(y.dat), 
                               y = as.numeric(y.dat$igr.mu),
                               t= as.numeric(y.dat$t),
                               R=as.numeric(y.dat$z1),
                               tem=as.numeric(y.dat$z2),
                               Npop =length(x.id.e),
                               Ndiet =length(x.id.food), #
                               Index=y.dat$Index.e,
                               z=y.dat$Index.food, #
                               obs_error_i=y.dat$igr.sd,
                               nsamp=sqrt(y.dat$nsamp),
                               # prior for mu and sigma when is shared
                               pr_mean_mu=log(20), pr_mean_mu_sd=log(1.5),
                               pr_mean_sigma= log(15), pr_mean_sigma_sd=log(10),
                               # prior for mu and sigma when hierarchical
                               pr_hyper_mu=20, pr_hyper_mu_sd=2.5,
                               pr_hyper_mu_deviation=1.5,pr_hyper_mu_deviation_sd=0.5,
                               pr_hyper_sigma=15,   pr_hyper_sigma_sd=10,
                               pr_hyper_sigma_deviation=1.5,pr_hyper_sigma_deviation_sd=0.5,
                               # priors for covariates
                               pr_a=log(1), pr_a_sd= 0.5, # updated
                               pr_bpow=log(1), pr_bpow_sd=0.5,
                               pr_bquad=log(1),pr_bquad_sd=0.5,
                               pr_diet=0, pr_diet_sd=1, 
                               # prior for observation error
                               pr_sigma_obs=5, pr_sigma_obs_sd=3)

### compile models: only the best models are compiled for space reason (each compiled file is > 6 MB)
#M1_0_J.mod <-  cmdstan_model(file.path(models.dir, 'M10J.stan'))
#M1_1_J.mod <-  cmdstan_model(file.path(models.dir, 'M11J.stan'))
#M1_2_J.mod <-  cmdstan_model(file.path(models.dir, 'M12J.stan'))
#M1_3_J.mod <-  cmdstan_model(file.path( models.dir, 'M13J.stan'))
M1_4_J.mod <-  cmdstan_model(file.path( models.dir, 'M14J.stan'))
#M1_0_I.mod <-  cmdstan_model(file.path( models.dir, 'M10I.stan'))
#M1_1_I.mod <-  cmdstan_model(file.path( models.dir, 'M11I.stan'))
#M1_2_I.mod <-  cmdstan_model(file.path( models.dir, 'M12I.stan'))
#M1_3_I.mod <-  cmdstan_model(file.path( models.dir, 'M13I.stan'))
#M1_4_I.mod <-  cmdstan_model(file.path( models.dir, 'M14I.stan'))
#M1_0_N.mod <-  cmdstan_model(file.path( models.dir, 'M10N.stan'))
#M1_1_N.mod <-  cmdstan_model(file.path( models.dir, 'M11N.stan'))
#M1_2_N.mod <-  cmdstan_model(file.path( models.dir, 'M12N.stan'))
#M1_3_N.mod <-  cmdstan_model(file.path( models.dir, 'M13N.stan'))
#M1_4_N.mod <-  cmdstan_model(file.path( models.dir, 'M14N.stan'))

### run models
x.elements=ls()
x.elements=x.elements[grep('M1',x.elements)]
x.elements=x.elements[-grep('data',x.elements)]
x.list=lapply(x.elements, get)
names(x.list)=x.elements
for(x in 1:length(x.list)){
  # get model type
  
  x.name=str_remove(names(x.list)[x], '.mod')
  x.type=substr(x.name,6,6)
  
  # create initial values
  if(x.type=='J'){
    x.init=list('sigma'=rep(20,length(x.id.e)+1),
                'mean_sigma'=20,
                'sigma_sigma'=6,
                'log_mu'=2,
                'sigma_obs'=3,
                'log_b1'=0,
                'log_b2'=1.2)
  }else if(x.type=='I'){
    x.init=list('mu'=rep(15,length(x.id.e)+1),
                'sigma_obs'=3,
                'log_b1'=0,
                'log_b2'=1.2)
  }else{
    x.init=NULL
  }
  # get model
  x.mod=x.list[[x]]
  # sample
  if(x.type =='N'){
    x.fit <- x.mod$sample(iter_warmup = 5000, 
                          iter_sampling = 5000,
                          data = data.M1,
                          seed = 123,
                          chains = 4,
                          parallel_chains = 4,
                          adapt_delta = 0.95,
                          max_treedepth=15,
                          thin=10,
                          show_exceptions = F)
  }else{
    x.fit <- x.mod$sample(iter_warmup = 5000, 
                          iter_sampling = 5000,
                          data = data.M1,
                          seed = 123,
                          chains = 4,
                          parallel_chains = 4,
                          adapt_delta = 0.95,
                          thin=10,
                          max_treedepth=15,
                          init = list(x.init,x.init,x.init,x.init),
                          show_exceptions = F)  
  }
  # save
  fit_mcmc=x.fit
  x.fit$save_output_files(dir = save.dir, 
                           basename = x.name, timestamp = FALSE, random = FALSE)
  saveRDS(x.fit, file.path(save.dir, paste0(x.name, '_fit.RDS')))
}


### M2 models ####
y.dat.a=x.dat[x.dat$t>16,]
set.seed(14)
y.dat.a=y.dat.a[sample(nrow(y.dat.a), nrow(y.dat.a)),]

# create indices for hierarchy ####
y.dat=y.dat.a
x.id.t=unique(y.dat$exp_id)
y.dat$Index.t=NA
for(i in 1:nrow(y.dat)){
  i.id=y.dat[i,]$exp_id
  y.dat[i,]$Index.t=which(x.id.t==i.id)
}

x.id.e=unique(paste0(y.dat$experiment))
y.dat$Index.e=NA
for(i in 1:nrow(y.dat)){
  i.id=paste0(y.dat[i,]$experiment)
  y.dat[i,]$Index.e=which(x.id.e==i.id)
}

x.id.no=1
y.dat$Index.no=1

# create indices for food items ####
x.id.food=unique(y.dat$main_prey)
y.dat$Index.food=NA
for(i in 1:nrow(y.dat)){
  i.id=y.dat[i,]$main_prey
  y.dat[i,]$Index.food=which(x.id.food==i.id)
}

write.csv(y.dat, file.path('data', 'M2_model_data.csv'), row.names = F)

# compile models: only the best models are compiled for space reason (each compiled file is > 6 MB)
#M2_0_Z.mod=cmdstan_model(file.path( models.dir, 'M20Z.stan'))
#M2_1_Z.mod=cmdstan_model(file.path( models.dir, 'M21Z.stan'))
#M2_2_Z.mod=cmdstan_model(file.path( models.dir, 'M22Z.stan'))
#M2_3_Z.mod=cmdstan_model(file.path( models.dir, 'M23Z.stan'))
M2_4_Z.mod=cmdstan_model(file.path( models.dir, 'M24Z.stan'))
#M2_0_K.mod=cmdstan_model(file.path( models.dir, 'M20K.stan'))
#M2_1_K.mod=cmdstan_model(file.path( models.dir, 'M21K.stan'))
#M2_2_K.mod=cmdstan_model(file.path( models.dir, 'M22K.stan'))
#M2_3_K.mod=cmdstan_model(file.path( models.dir, 'M23K.stan'))
#M2_4_K.mod=cmdstan_model(file.path( models.dir, 'M24K.stan'))
#M2_0_N.mod=cmdstan_model(file.path( models.dir, 'M20N.stan'))
#M2_1_N.mod=cmdstan_model(file.path( models.dir, 'M21N.stan'))
#M2_2_N.mod=cmdstan_model(file.path( models.dir, 'M22N.stan'))
#M2_3_N.mod=cmdstan_model(file.path( models.dir, 'M23N.stan'))
#M2_4_N.mod=cmdstan_model(file.path( models.dir, 'M24N.stan'))

data.M2 <- list(N = nrow(y.dat), 
                                     y = as.numeric(y.dat$igr.mu),
                                     t= as.numeric(y.dat$t),
                                     R=as.numeric(y.dat$z1),
                                     tem=as.numeric(y.dat$z2),
                                     Npop =length(x.id.e),
                                     Ndiet =length(x.id.food), #
                                     Index=y.dat$Index.e,
                                     z=y.dat$Index.food, #
                                     obs_error_i=y.dat$igr.sd,
                                     nsamp=sqrt(y.dat$nsamp),
                                     # prior for mu and sigma when is shared
                                     pr_mean_A=10, pr_mean_A_sd=10,
                                     pr_mean_B= 10, pr_mean_B_sd=10,
                                     # prior for mu and sigma when hierarchical
                                     pr_hyper_A=5, pr_hyper_A_sd=1.5,
                                     pr_hyper_A_deviation=1.5,pr_hyper_A_deviation_sd=1,
                                     pr_hyper_B=10,   pr_hyper_B_sd=1.5,
                                     pr_hyper_B_deviation=1.5,pr_hyper_B_deviation_sd=1,
                                     # priors for covariates
                                     pr_C=log(0.05), pr_C_sd= 0.5,
                                     pr_bpow=log(1), pr_bpow_sd=log(2),
                                     pr_bquad=log(3),pr_bquad_sd=log(2),
                                     pr_diet=log(0.5), pr_diet_sd=0.5, 
                                     # prior for observation error
                                     pr_sigma_obs=5, pr_sigma_obs_sd=3)

x.elements=ls()
x.elements=x.elements[grep('M2',x.elements)]
x.elements=x.elements[-grep('data',x.elements)]
x.list=lapply(x.elements, get)
names(x.list)=x.elements
for(x in 1:length(x.list)){
  # get model type
  x.name=str_remove(names(x.list)[x], '.mod')
  x.type=substr(x.name,6,6)
  # get model
  x.mod=x.list[[x]]
  # sample
  x.fit <- x.mod$sample(iter_warmup = 5000, 
                        iter_sampling = 5000,
                        data = data.M2,
                        seed = 123,
                        chains = 4,
                        parallel_chains = 4,
                        adapt_delta = 0.95,
                        max_treedepth=15,
                        thin=10,
                        show_exceptions = F)
# save
  x.fit$save_output_files(dir = save.dir, 
                        basename = x.name, timestamp = FALSE, random = FALSE)
  saveRDS(x.fit, file.path(save.dir, paste0(x.name, '_fit.RDS')))
}


