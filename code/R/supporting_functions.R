### write here some description

# growth models ####
un.gr.formula=function(w1,gr,t1,t2){exp(log(w1)+((t2-t1)/100*(gr)))}

gr.formula=function(w1,w2,t1,t2){((log(w2)-log(w1))*100)/(t2-t1)}

# general functions ####
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

lnorm.transformation=function(mu.w, sd.w, out='mn', npop=500){
  w.sigma.hat=sqrt(log((sd.w^2 / (mu.w^2)) + 1))
  w.mu.hat=log(mu.w) - (w.sigma.hat^2) / 2
  if(out=='distr'){
    w.obs=rlnorm(npop, meanlog = w.mu.hat, sdlog = w.sigma.hat)    
    return(w.obs)
  }else{
    return(c(w.mu.hat, w.sigma.hat))
  }
}

growth.pop.stat=function(w1,w2,t1,t2){
  # Compute the covariance matrix between log w1 and w2
  cov_matrix <- cov(cbind(log(w1), log(w2)))
  # Extract variances and covariance
  var_log_w1 <- cov_matrix[1, 1]
  var_log_w2 <- cov_matrix[2, 2]
  cov_log_w1_w2 <- cov_matrix[1, 2]
  # Growth rate formula
  igr.mu <- (log(mean(w2)) - log(mean(w1))) * 100 / (t2 - t1)
  # Standard deviation of the growth rate using covariance adjustment
  igr.sd <- sqrt((var_log_w2 + var_log_w1 - 2 * cov_log_w1_w2) * (100^2) / (t2 - t1)^2)
  return(c(igr.mu, igr.sd))
}

sigmoid_transition <- function(t, m1,m2, x1,x2, k, type='sigm') {
  if(type=='sigm'){
    # Define Gaussian and Exponential functions
    # Sigmoid transition function
    phi0=data.frame(t=t,phi=NA)
    phi0$phi=ifelse(phi0$t<x1,0,
                    ifelse(phi0$t>x2,1,
                           1 / (1 + exp(-0.25 * (t - (x1+x2)/2)))))
    phi=phi0$phi
    # Combine functions with sigmoid transition
    (1 - phi) * m1 + phi * m2 
    
  }else{
    phi <- pmax(pmin((t - (x - k)) / (2 * k), 1), 0)
    # Combine functions with smooth transition
    ifelse(t < x - k, mj,
           ifelse(t > x + k, m2,
                  (1 - phi) * m1 + phi * m2))
  }
}

apply_sigmoid_row <- function(row1, row2, t_values) {
  mapply(function(m1, m2, t) {
    sigmoid_transition(t = t, m1 = m1, m2 = m2, x1 = 16, x2 = 41, k = 0.2)
  }, row1, row2, t_values)
}


format.posteriors=function(x.model){
  x.post=x.model$draws()
  x.post=as.data.frame(x.post)
  x.post=x.post[,-grep('y.rep', names(x.post))]
  x.post=x.post[,-grep('raw', names(x.post))]
  x.post=x.post[,-grep('log_lik', names(x.post))]
  x.post=x.post[,-grep('lp_', names(x.post))]
  x.vars=unique(substr(names(x.post),3, nchar(names(x.post))))
  n.chain=ncol(x.post)/length(x.vars)
  x.post.matrix=matrix(ncol = length(x.vars), nrow=nrow(x.post)*n.chain)
  for(i in 1:length(x.vars)){
    x.post.matrix[,i]=as.numeric(unlist(x.post[,substr(colnames(x.post),3, nchar(colnames(x.post)))==x.vars[i]]  ))
  }
  x.post.matrix=as.data.frame(x.post.matrix)
  names(x.post.matrix)=x.vars
  return(x.post.matrix)
}


## best models formulas ####

M1=function(a,b1,b2,b3,sigma,mu,t,tem,R,Fi){
  T_t=tem+0.5
  resp= (a * (R*Fi)^b1* (b2*T_t+b3*T_t^2)) * exp(-(t - mu)^2 / ( 2 * sigma^2))
  return(resp)
}
M2=function(A,b1,b4a,b2,b3,B,C,t,tem,R,Fi){
  T_t=tem+0.5
  resp=A*(T_t^b4a)*(R*Fi)^b1 + (B*(b2*T_t + b3*T_t^2) * (R*Fi)^b1) * exp(-t * C)
  return(resp)
}

# estimation models ####
# EM1
EM1=function(w.matrix, biol.value){
  w.matrix$w.mu.pop.transf=exp(lnorm.transformation(mu.w =w.matrix$w.mu.pop, sd.w = w.matrix$w.sd.pop )[1:(nrow(w.matrix))])
  
  output=data.frame(gr.pred=gr.formula(w1=w.matrix[1:(nrow(biol.value)),]$w.mu.pop.transf,
                                       w2=w.matrix[2:(nrow(biol.value)+1),]$w.mu.pop.transf,
                                       t1=biol.value$t1,
                                       t2=biol.value$t2),
                    gr.obs=biol.value$igr.mu,
                    t=biol.value$t1)
  
  
  igr.sd.stat=NULL
  for(i in 1:(nrow(w.matrix)-1)){
    var1=log(1+(w.matrix[i,]$w.sd.pop^2)/(w.matrix[i,]$w.mu.pop^2))
    var2=log(1+(w.matrix[i+1,]$w.sd.pop^2)/(w.matrix[i+1,]$w.mu.pop^2))
    pred.sd=sqrt(abs(var2-var1)*((100^2)/((biol.value[i,]$t2-biol.value[i,]$t1)^2)))
    igr.sd.stat=rbind(igr.sd.stat, pred.sd)
  }  
  output$igr.sd.obs=biol.value$igr.sd
  output$igr.sd.pred=igr.sd.stat[,1]
  return(output)
}

# EM2
EM2=function(w.matrix, biol.value, s.size, iter){
  sample.size=s.size
  biol.pars
  k.res=NULL
  for(k in 1:(nrow(w.matrix)-1)){
    mu.w1=w.matrix[k,]$w.mu.pop
    sd.w1=w.matrix[k,]$w.sd.pop
    mu.w2=w.matrix[k+1,]$w.mu.pop
    sd.w2=w.matrix[k+1,]$w.sd.pop
    
    # transform values
    w1.sigma.hat=sqrt(log((sd.w1^2 / (mu.w1^2)) + 1))
    w1.mu.hat=log(mu.w1) - (w1.sigma.hat^2) / 2
    w2.sigma.hat=sqrt(log((sd.w2^2 / (mu.w2^2)) + 1))
    w2.mu.hat=log(mu.w2) - (w2.sigma.hat^2) / 2
    
    
    res=NULL
    for(j in 1:iter){
      w1.obs=rlnorm(s.size, meanlog = w1.mu.hat, sdlog = w1.sigma.hat)
      w2.obs=rlnorm(s.size, meanlog = w2.mu.hat, sdlog = w2.sigma.hat)
      
      gr.samples=(gr.formula(w2=w2.obs,w1=w1.obs,t2=biol.value[k,]$t2, t1=biol.value[k,]$t1))
      gr.samples.s=(gr.formula(w2=sort(w2.obs),w1=sort(w1.obs),t2=biol.value[k,]$t2, t1=biol.value[k,]$t1))
      
      i.res=data.frame(mu=mean(gr.samples),sd=sd(gr.samples),
                       mu.s=mean(gr.samples.s),sd.s=sd(gr.samples.s))
      res=rbind(res, i.res)
    }
    res=res%>%
      dplyr::summarise(mu.mu=mean(mu), mu.sd=sd(mu), sd.mu=mean(sd), sd.sd=sd(sd),
                       mu.mu.s=mean(mu.s), mu.sd.s=sd(mu.s), sd.mu.s=mean(sd.s), sd.sd.s=sd(sd.s),.groups = "keep")
    k.res=rbind(k.res, res)
  }
  
  return(k.res)
  
}

# EM3
EM3=function(w.matrix, biol.value, s.size, iter, sd.vector){
  igr.sd.sens=NULL
  for(i in 1:(nrow(w.matrix)-1)){
    mu.w1=w.matrix[i,]$w.mu.pop
    sd.w1=w.matrix[i,]$w.sd.pop
    mu.w2=w.matrix[i+1,]$w.mu.pop
    sd.w2=w.matrix[i+1,]$w.sd.pop
    
    # transform values
    w1.sigma.hat=sqrt(log((sd.w1^2 / (mu.w1^2)) + 1))
    w1.mu.hat=log(mu.w1) - (w1.sigma.hat^2) / 2
    w2.sigma.hat=sqrt(log((sd.w2^2 / (mu.w2^2)) + 1))
    w2.mu.hat=log(mu.w2) - (w2.sigma.hat^2) / 2
    
    # test vector
    sd.matrix=data.frame(sd.val=sd.vector, p.val=NA)
    for(j in 1:length(sd.vector)){
      p.vec=NULL
      for(k in 1:iter){
        # generate
        w1.obs=rlnorm(1000, meanlog = w1.mu.hat, sdlog = w1.sigma.hat)
        w2.obs=rlnorm(1000, meanlog = w2.mu.hat, sdlog = w2.sigma.hat)
        gr.obs=rnorm(n=1000, mean=biol.value[i,]$igr.mu, sd=sd.vector[j])
        w2.pred=(un.gr.formula(w1=w1.obs, gr=gr.obs, t1=biol.value[i,]$t1, t2=biol.value[i,]$t2))
        j.test=ks.test((w2.pred), (w2.obs))
        p.vec=c(p.vec,j.test$p.value )
      }
      sd.matrix[j,]$p.val=mean(p.vec)
    }
    pred.sens=sd.matrix[sd.matrix$p.val==max(sd.matrix$p.val),]$sd.val
    igr.sd.sens=rbind(igr.sd.sens, pred.sens)
  }
  return(igr.sd.sens)
}

# compare distributions
# Equations S7-S10
testing.predictions=function(xdat, xobs, xvar, alpha=0.05){
  times.vec=unique(xdat$t)
  significance=matrix(nrow = length(times.vec))
  for(tt in 1:length(times.vec)){
    Xj=xdat[xdat$t==times.vec[tt],xobs]
    Yj=xdat[xdat$t==times.vec[tt],xvar]
    Wj=Xj-Yj
    Xhat=mean(Xj)
    Yhat=mean(Yj)
    What=Xhat-Yhat
    VARW=sum((Wj-What)^2)/(length(Xj)*(length(Yj)-1))
    sign.int=qt(1 - alpha/2, df = length(Yj) - 1)
    VAR_int=sqrt(VARW)*sign.int
    ci.lo=What-VAR_int
    ci.hi=What+VAR_int
    if(sign(ci.lo)==sign(ci.hi)){
      significance[tt,1]=1    
    }else{
      significance[tt,1]=0
    }
  }
  return(significance)
}

#### some common plots ####
plot.sample.sd=function(xdat, x.title, x.color, x.var){
  xdat=cbind(xdat[, c("t", "real.mu", "real.sd",'obs')],
        xdat[,grep(x.var, names(xdat))])
  names(xdat)[5:6]=c('pred.mu','pred.sd')
 p=ggplot(data=xdat)+
   ylab(expression(sigma[IGR][t]^2))+
  xlab('Time (days)')+
    
    geom_point(aes(x=t+0.5, y=pred.mu), color=x.color)+
    geom_errorbar(aes(x=t+0.5, ymin=pred.mu-pred.sd,
                      ymax=pred.mu+pred.sd), color=x.color)+
    geom_point(aes(x=t, y=obs),shape=18)+
    geom_point(aes(x=t, y=real.mu),shape=18, color='gold')+
    geom_errorbar(aes(x=t, ymin=real.mu-real.sd,
                      ymax=real.mu+real.sd), color='black')+
    
    ggtitle(x.title)
  return(p)
}

base.plot.curve=function(xdat){
  ggplot(data=xdat, aes(x=xvar, y=resp.mu))+
    geom_line()+
    geom_ribbon(aes(ymin=resp.lo, ymax=resp.hi), alpha=0.2)+
    theme_minimal()
}



### models diagnostics ####

# prior posterior
posterior.diags=function(fit_mcmc, prior.list, fixed.pars, label_list, random.pars,label,save.dir){
  posterior_draws = fit_mcmc$draws()
  post.dr=as.data.frame(posterior_draws)
  x.post=as.data.frame(posterior_draws)
  x.post=x.post[,-grep('y_rep', names(x.post))]
  x.post=x.post[,-grep('lik', names(x.post))]
  sel.vector=seq(1, ncol(x.post),4)
  numbs=seq(0,3,1)
  w.matrix=matrix(ncol=ncol(x.post)/4, nrow=2000, NA)
  w.matrix=as.data.frame(w.matrix)
  for(i in 1:length(sel.vector)){
    base.sel=sel.vector[i]
    i.sel=base.sel+numbs
    i.values=x.post[,i.sel]
    i.nm=unique(substr(names(x.post[,i.sel]),3, unique(nchar(names(x.post[,i.sel])))))
    i.values=unlist(i.values,use.names = FALSE)
    w.matrix[,i]=i.values
    names(w.matrix)[i]=i.nm
  }
  pars.list=data.frame(par=c(fixed.pars, random.pars[c(1:2)] ))
  pars.list$type='fix'
  pars.list[(nrow(pars.list)-1):nrow(pars.list),]$type='random'
  pars.list=pars.list[pars.list$par%in% names(w.matrix),]
  plot.list=list()
  for(x in 1:nrow(pars.list)){
    x.par=pars.list[x,]
    if(x.par$type=='fix'){
      if(x.par$par %in% c('sigma_obs', 'A','B', 'bTA')){
        x.post=w.matrix[,names(w.matrix)==paste0( x.par$par)]
      }else{
        x.post=w.matrix[,names(w.matrix)==paste0(x.par$par)]
      }
      if(x.par$par%in%c('A','B')){x.par$par=paste0('mean_',x.par$par)}
      if(x.par$par%in%c('mu','sigma')){x.par$par=paste0('mean_',x.par$par)}
      
    }else if(x.par$type=='random'){
      x.post=w.matrix[,names(w.matrix)==paste0( x.par$par)]
      if(substr(x.par$par,1,3)=='mea'){x.par$par=str_replace(x.par$par, 'mean', 'hyper')}
      if(substr(x.par$par,1,3)=='sig'){
        if(substr(x.par$par,7,8)=='si'){x.par$par='hyper_sigma_deviation'}
        if(substr(x.par$par,7,8)=='mu'){x.par$par='hyper_mu_deviation'}
        if(substr(x.par$par,7,7)=='A'){x.par$par='hyper_A_deviation'}
        if(substr(x.par$par,7,7)=='B'){x.par$par='hyper_B_deviation'}
      }
    }
    x.post=data.frame(type='posterior', value=x.post)
    if(x.par$par %in% c('bF','bT','bTA','bTB')){
      prior.pars=prior.list[,grep(paste0('pr_','bpo'), names(prior.list))]
    }else if(x.par$par %in% c('bT2','bT1','bTB1', 'bTB2')){
      prior.pars=prior.list[,grep(paste0('pr_','bqu'), names(prior.list))]
    }else if(x.par$par %in% c('a1','a2')){
      prior.pars=prior.list[,grep(paste0('pr_','a'), names(prior.list))]
    }else{
      prior.pars=prior.list[,grep(paste0('pr_', x.par$par), names(prior.list))]
    }
    x.prior=data.frame(type='prior', value= rnorm(n=2000, mean=prior.pars[[1]], sd=prior.pars[[2]]))
    if(x.par$par=='bTB2'){
      x.prior=data.frame(type='prior', value= rnorm(n=2000, mean=-prior.pars[[1]], sd=prior.pars[[2]]))
    }
    i.dist=rbind(x.prior, x.post)
    i.pl=ggplot(data=i.dist, aes(value))+
      geom_density(aes(fill=type), alpha=0.4)+
      ggtitle(label_list[[x]])+
      scale_fill_viridis_d()+
      xlab('Value')+
      ylab('Density')+
      theme(legend.position = 'bottom')+
      labs(fill='')
    plot.list[[x]]=i.pl
  }
  pl=ggpubr::ggarrange(plotlist = plot.list, common.legend = T, legend = 'bottom')
  ggsave(plot=pl,file.path(save.dir, model.lab, paste0(model.lab,'_ppdist.png')), width = 20, height = 20, units='cm')
}

# mcmc diags ###
mcmc.diags=function(x.fit, x.post, np_cp, lp_cp, label, fixed.pars, random.pars, n.iter=2000, version=NULL){
  image.size=length(fixed.pars)
  if(image.size<=3){
    x.wid=20
    x.hei=8
  }else if(image.size %in% 4:6){
    x.wid=20
    x.hei=16
  }else if(image.size %in% 6:12){
    x.wid=25
    x.hei=24
  }else if(image.size > 12){
    x.wid=30
    x.hei=32
  }
  
  
  if(version=='paper'){
    
 
  iterations <- 200:500
  breaks_every_10 <- iterations[seq(1, length(iterations), by = 150)]
  param_names=c(fixed.pars,random.pars)
  param_labels <- setNames(
    str_replace_all(param_names, c(
      "sigma_obs" = "sigma[obs]",   
      "sigma_sigma" = "Hyper~sigma[sigma]",
      "mean_sigma" = "Hyper~sigma[mu]",
      "sigma_B" = "Hyper~B[sigma]",
      "mean_B" = "Hyper~B[mu]",
      "bF" = "b[F]",
      "bT1" = "b[T1]",   
      "bTA" = "b[TA]", 
      "bTB1" = "b[TB1]", 
      "bTB2" = "b[TB2]", 
      "bT2" = "b[T2]"  
    )),
    param_names
  )
  param_labels <- setNames(
    str_replace(param_labels, "sigma_(\\d+)", "sigma[\\1]"),
    param_names
  )
  param_labels <- as_labeller(param_labels, label_parsed)
 
  p=mcmc_trace(x.post, pars =  c(fixed.pars,random.pars), n_warmup = 250, window = c(200, 500))+
    facet_wrap(~parameter, labeller = param_labels, scales = 'free')
  p=p+  scale_x_continuous(breaks = unique(breaks_every_10), labels = 10*(unique(breaks_every_10)))
  p
  png(file.path(save.dir, label, paste0(label,'_trace_random.png')), width = 25, height = 27, units='cm', res=300)
  p+theme_minimal()
  dev.off()
  }
 
  param_labels <- setNames(
    str_replace_all(param_names, c(
      "sigma_obs" = "sigma[obs]",   
      "sigma_sigma" = "Hyper~sigma[sigma]",
      "mean_sigma" = "Hyper~sigma[mu]",
      "sigma_B" = "Hyper~B[sigma]",
      "mean_B" = "Hyper~B[mu]",
      "bF" = "b[F]",
      "bT1" = "b[T1]",   
      "bTA" = "b[TA]", 
      "bTB1" = "b[TB1]", 
      "bTB2" = "b[TB2]", 
      "bT2" = "b[T2]"  
    )),
    param_names
  )
  
  
  # get stats
  ess=neff_ratio(x.fit)
  rhats=rhat(x.fit)
  x.diags=x.fit$sampler_diagnostics()
  x.diags=x.diags[,,2]
  x.diags=as.data.frame(x.diags)
  div.iter=sum(apply(x.diags, 1, sum))
  div.iter.perc=div.iter/n.iter
  mcmc.stats=data.frame(rhat.max=max(rhats), ess.min=min(ess), no.divg=(div.iter.perc))
  return(mcmc.stats)
}




rmse.fit=function(posterior_draws,y.dat,save.dir,label,x.fit, niter=2000){
  x.post=as.data.frame(posterior_draws)
  w.pred=x.post[,grep('y_rep', names(x.post))]
  sel.vector=seq(1, ncol(w.pred),4)
  numbs=seq(0,3,1)
  w.matrix=matrix(ncol=nrow(y.dat), nrow=niter, NA)
  for(i in 1:length(sel.vector)){
    base.sel=sel.vector[i]
    i.sel=base.sel+numbs
    i.values=w.pred[,i.sel]
    i.values=unlist(i.values,use.names = FALSE)
    w.matrix[,i]=i.values
  }
  w.pred=w.matrix
  i.dat=y.dat
  matrix.0=matrix(nrow=nrow(w.pred), ncol=ncol(w.pred))
  for(i in 1:nrow(i.dat)){
    # posteriors have obs by column and iteration by row
    pi=w.pred[,i]
    oi=i.dat[i,]$igr.mu
    matrix.0[,i]=pi-oi
  }
  rmse.store=vector(length = nrow(matrix.0))# rmse distribution
  for(i in 1:nrow(matrix.0)){
    # posteriors have obs by column and iteration by row
    rmse.store[i]=sqrt(sum(matrix.0[i,]^2)/ncol(matrix.0))
  }
  quants <- matrix(0,nrow=3,ncol=ncol(w.pred))# predictions
  for (i in 1:ncol(w.pred))
    quants[,i] <- quantile(w.pred[,i],probs=,c(0.05,0.5,0.95))
  
  i.dat$pred=quants[2,]
  i.dat$lo.pred=quants[1,]
  i.dat$hi.pred=quants[3,]
  i.dat$p.err=i.dat$pred-i.dat$igr.mu
  
  # residuals
  baseplot=ggplot(data=i.dat)+
    ylab('Residuals')
  
  p1=baseplot+
    geom_point(aes(x=t, y=p.err))+
    xlab('Time (DAH)')+
    geom_hline(yintercept = 0, color='red')
  p2=baseplot+
    geom_point(aes(x=Temperature, y=p.err))+
    xlab('Temperature')+
    geom_hline(yintercept = 0, color='red')
  p3=baseplot+
    geom_point(aes(x=Rcomma, y=p.err))+
    xlab("R'")+
    geom_hline(yintercept = 0, color='red')
  p4=baseplot+
    geom_point(aes(x=sample(1:nrow(y.dat),nrow(y.dat)), y=p.err))+
    xlab('Id')+
    geom_hline(yintercept = 0, color='red')
  p5=ggpubr::ggarrange(p1,p2,p3,p4)
  ggsave(plot=p5,paste0(save.dir, label,'_residuals.png'), width = 20, height = 15, units='cm')
  
  # fits
  y.pred=i.dat%>%
    dplyr::group_by(t,Temperature,Rcomma,exp_id)%>%
    dplyr::summarise(.groups="keep",pred=median(pred), lo.pred=median(lo.pred),hi.pred=median(hi.pred))
  
  # fit to data
  summary.data=data.frame(table(y.pred$exp_id, y.pred$t))
  summary.data=summary.data[summary.data$Freq>0,]
  summary.data=data.frame(table(summary.data$Var1))
  summary.data=summary.data[summary.data$Freq==1,]
  dummy.add=y.pred[y.pred$exp_id %in% summary.data$Var1,]
  dummy.add$t=dummy.add$t+10
  dummy.less=y.pred[y.pred$exp_id %in% summary.data$Var1,]
  dummy.less$t=dummy.less$t-10
  y.pred=rbind(y.pred, dummy.add, dummy.less)
  y.dat$Temperature=paste0(y.dat$Temperature, ' °C')
  y.pred$Temperature=paste0(y.pred$Temperature, ' °C')
  
  p.fit=ggplot(data=y.dat, aes(group=1))+
    facet_wrap(~Temperature+exp_id, ncol=10)+
    geom_point(aes(x=t, y=igr.mu))+
    geom_line(data=y.pred,aes(x=t, y=pred, group = 1), color='red')+
    theme(legend.position = 'bottom')+
    #ggtitle('Red line is model fit')+
    scale_fill_viridis_c()+
    xlab('Time (DAH)')+
    labs(fill= "R'")+
    geom_ribbon(data=y.pred,aes(x=t, ymin=lo.pred, ymax=hi.pred, fill=Rcomma, group = 1), alpha=0.3)
  ggsave(plot=p.fit,paste0(save.dir, label,'_fit.png'), width = 30, height = 30, units='cm')
  
  log_lik <- x.fit$draws(variables = "log_lik")
  loo_result=x.fit$loo()
  
  png(paste0(save.dir, label,'_paretok.png'), width = 15, height = 15, units='cm', res=500)
  print(plot(loo_result))
  dev.off()
  pareto_k=loo_result$diagnostics$pareto_k
  pareto_k=length(pareto_k[pareto_k>=0.7])
  loo_ests=as.data.frame(loo_result$estimate)
  colnames(loo_ests)=tolower(colnames(loo_ests))
  loo_ests=rbind(loo_ests,data.frame(estimate=pareto_k, se=NA, row.names = 'pareto.k'))
  loo_ests=rbind(loo_ests,data.frame(estimate=mean(rmse.store), se=sd(rmse.store), row.names = 'rmse'))
  loo_ests$par=tolower(rownames(loo_ests))
  return(loo_ests)
}


summarise.metrics=function(store,save.dir, save.lab){
  
  xstore=store%>%replace(is.na(.),0)
  xstore=xstore[xstore$par %in% c( "ess.min" ,  "no.divg"  , "rhat.max", "pareto.k","rmse","elpd_loo"),]
  legend=data.frame(par=unique(xstore$par),type=c('mcmc','mcmc','mcmc','pred.power','pred.power','pred.power'))
  xstore=left_join(xstore, legend, by='par')
  xstore$model.id=substr(xstore$model,4,4)
  #xstore$model.id=str_remove(xstore$model.id, "\\.$"  )
  xstore$model.h=substr(xstore$model,6,nchar(xstore$model))
  labs=c("ess.min" ='a) Lowest Ess',  "no.divg" ='b) % of Divergent transitions' ,"rhat.max"= "c) Largest Rhat")
  p1=ggplot(data=xstore[xstore$type=='mcmc',], aes(x=reorder(model.id, as.numeric(str_remove(model.id,'m'))), y=estimate, fill=model.h))+
    geom_col(position='dodge')+
    #scale_y_continuous(expand=c(0,0))+
    #geom_errorbar(aes(ymin=(estimate-se),ymax=(estimate+se)))+
    coord_cartesian(ylim = c(min(xstore[xstore$type=='mcmc',]$estimate), max(xstore[xstore$type=='mcmc',]$estimate)))+
    facet_wrap(~par, scales='free_y', labeller = labeller(par = labs))+
    #ggtitle('MCMC')+
    xlab('Fixed effect')+
    ylab('Value')+
    theme_minimal()+
    
    labs(color='Random effect')+
    scale_fill_manual(values=c("#31688EFF","#440154FF", "#6DCD59FF" ));p1
  labs=c( "pareto.k" = 'e) Pareto K',"rmse" = 'f) RMSE',"elpd_loo" ='d) elpd loo')
  p2=ggplot(data=xstore[xstore$type=='pred.power' ,], aes(x=reorder(model.id, as.numeric(str_remove(model.id,'m'))), y=estimate, color=model.h))+
    geom_point()+
    geom_errorbar(aes(ymin=(estimate-(se)),ymax=(estimate+(se))), 
                  width=0.15)+
    facet_wrap(~par, scales='free', labeller = labeller(par = labs))+
    #ggtitle('Predictive Accuracy')+
    theme_minimal()+
    xlab('Fixed effect')+
    ylab('Value')+
    labs(color='Random effect')+
    scale_color_manual(values=c("#31688EFF","#440154FF", "#6DCD59FF" )); p2
  
  #xstore[xstore$par=='p_loo',]$estimate=xstore[xstore$par=='p_loo',]$estimate/xstore[xstore$par=='p_loo',]$n.parameters
  ggpubr::ggarrange(p1,p2,ncol=1, common.legend = T)
  ggsave(file.path(save.dir, paste0(save.lab,'.jpeg')), width = 30, height = 15, units='cm')
  
}














