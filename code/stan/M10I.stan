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
  real mean_mu;              // peak position (shift)
  real<lower=0> sigma_mu;  
  vector [Npop+1] raw_mu;      // Mean of the Gaussian curve for each population   
  real log_sigma; 
  real log_a;      // scale
  real sigma_obs;
}
 transformed parameters {
  vector[N] obs_error = (sigma_obs + obs_error_i)./nsamp; 
	real a=exp(log_a);      
	vector<lower=0> [Npop+1] mu=(mean_mu + sigma_mu * raw_mu);  
	real sigma = exp(log_sigma);   
}
model {
	//Intercept
	log_a ~ normal(pr_a, pr_a_sd);
	
	// Gaussian
	mean_mu ~ normal(pr_hyper_mu, pr_hyper_mu_sd);
	sigma_mu ~ normal(pr_hyper_mu_deviation, pr_hyper_mu_deviation_sd);
	raw_mu ~ std_normal();
	log_sigma ~ normal(pr_mean_sigma, pr_mean_sigma_sd);
	
	// Obs_error
	sigma_obs ~ normal(pr_sigma_obs,pr_sigma_obs_sd);
	
	// Fit 	
	for(i in 1:N)
	{
		y[i] ~ normal(a * exp(-(t[i] - mu[Index[i]])^2 / ( 2* sigma^2)), obs_error[i]);
	}	
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  // Log likelihood for model comparisons
  for (i in 1:N){
    log_lik[i] = normal_lpdf(y[i] | a * exp(-(t[i] - mu[Index[i]])^2 / ( 2* sigma^2)), obs_error[i]);
 }
  // Generate posterior predictive checks
  for (i in 1:N) {
y_rep[i] = normal_rng( a * exp(-(t[i] - mu[Index[i]])^2 / ( 2* sigma^2)), obs_error[i]);
  }
}
