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
parameters {
  real mean_sigma;              // peak position (shift)
  vector [Npop+1] raw_sigma;      
  real<lower=0> sigma_sigma;  
  real log_mu; 
  real log_a;      // scale
  real sigma_obs;
}
 transformed parameters {
  vector[N] obs_error = (sigma_obs + obs_error_i)./nsamp; 
  real a=exp(log_a);      // scale
  vector<lower=0> [Npop+1] sigma=mean_sigma + sigma_sigma * raw_sigma;     
  real mu=exp(log_mu);  
}
model {
	log_a ~ normal(pr_a, pr_a_sd);
    
	mean_sigma ~ normal(pr_hyper_sigma, pr_hyper_sigma_sd);
	sigma_sigma ~ normal(pr_hyper_sigma_deviation, pr_hyper_sigma_deviation_sd);
	raw_sigma ~ std_normal();
	log_mu ~ normal(pr_mean_mu,pr_mean_mu_sd);
	
	sigma_obs ~ normal(pr_sigma_obs,5);
  	
	for(i in 1:N)
	{
		y[i] ~ normal( a * exp(-(t[i] - mu)^2 / ( 2* sigma[Index[i]]^2)), obs_error[i]);
	}	
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  // Log likelihood for model comparisons
  for (i in 1:N){
    log_lik[i] = normal_lpdf(y[i] |  a * exp(-(t[i] - mu)^2 / ( 2* sigma[Index[i]]^2)), obs_error[i]);
 }
  // Generate posterior predictive checks
  for (i in 1:N) {
y_rep[i] = normal_rng(  a * exp(-(t[i] - mu)^2 / ( 2* sigma[Index[i]]^2)), obs_error[i]);
  }
}
