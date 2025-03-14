data {
  // data
  int<lower=0> N;      
  int<lower=0> Npop; 
  array[N] int Index;
  int<lower=0> Ndiet; 
  array[N] int z; 
  vector[N] t;         
  vector[N] y;        
  vector[N] R;  	
  vector[N] tem;  
  vector[N] obs_error_i;
  vector[N] nsamp; 
  // priors
  real pr_a;
  real pr_a_sd;
  real pr_bpow;
  real pr_bpow_sd;
  real pr_bquad;
  real pr_bquad_sd;
  real pr_hyper_sigma;
  real pr_hyper_sigma_sd;
  real pr_hyper_sigma_deviation;
  real pr_hyper_sigma_deviation_sd;
  real pr_mean_sigma;
  real pr_mean_sigma_sd;
  real pr_hyper_mu;
  real pr_hyper_mu_sd;
  real pr_mean_mu;
  real pr_mean_mu_sd;
  real pr_sigma_obs;
  real pr_sigma_obs_sd;
  real pr_hyper_mu_deviation;
  real pr_hyper_mu_deviation_sd;
  real pr_diet;
  real pr_diet_sd;
}
transformed data {
  vector[N] R_t ;
  vector[N] T_t ;
	R_t = R+0.5;
	T_t = tem+0.5;
}
parameters {
  real log_mu;       // Mean of the Gaussian curve for each population
  real log_sigma;  // peak amplitude
  real log_a;      
  real sigma_obs;
  real log_bF;   
  real log_bT1;   
  real log_bT2;  
  vector [Ndiet-1] F_raw;     
}
 transformed parameters {
  vector[N] obs_error = (sigma_obs + obs_error_i)./nsamp; 
  real mu=exp(log_mu);      // scale
  real sigma = exp(log_sigma);
  vector [Ndiet] F = inv_logit(append_row(F_raw, -sum(F_raw)));    
  real a=exp(log_a);     
  real bF=exp(log_bF);   
  real bT1=exp(log_bT1);    
  real bT2=-exp(log_bT2);   
}
model {
  
	log_a ~ normal(pr_a, pr_a_sd);
	log_bF ~ normal(pr_bpow,pr_bpow_sd);
	log_bT1 ~ normal(pr_bquad, pr_bquad_sd);
	log_bT2 ~ normal(pr_bquad, pr_bquad_sd);
	
	log_sigma ~ normal(pr_mean_sigma, pr_mean_sigma_sd);
	log_mu ~ normal(pr_mean_mu,pr_mean_mu_sd);
	
	F_raw ~ normal(pr_diet, pr_diet_sd);
	sigma_obs ~ normal(pr_sigma_obs,pr_sigma_obs_sd);
	
	for(i in 1:N)
	{
		y[i] ~ normal( (a*(R[i]*F[z[i]])^bF * (bT1*T_t[i]+bT2*T_t[i]^2)) * exp(-(t[i] - mu)^2 / ( 2 * sigma^2)), obs_error[i]);
	}
}
generated quantities {
  vector[N] y_rep;
   vector[N] log_lik;
   // Log likelihood for model comparisons
   for (i in 1:N){
     log_lik[i] = normal_lpdf(y[i] |(a*(R[i]*F[z[i]])^bF * (bT1*T_t[i]+bT2*T_t[i]^2)) * exp(-(t[i] - mu)^2 / ( 2 * sigma^2)), obs_error[i]);
  } 
   // Generate posterior predictive checks
   for (i in 1:N) {
 y_rep[i] = normal_rng( (a*(R[i]*F[z[i]])^bF * (bT1*T_t[i]+bT2*T_t[i]^2)) * exp(-(t[i] - mu)^2 / ( 2 * sigma^2)), obs_error[i]);
   }
 
 }
