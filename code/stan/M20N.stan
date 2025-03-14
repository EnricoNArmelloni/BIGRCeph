data {
  int<lower=0> N;       // Number of observations
  int<lower=0> Npop; // number of populations
  array[N] int Index;
  array[N] int z; 
  int<lower=0> Ndiet; // number of diets
  vector[N] t;  
  vector[N] y;          // Response variable
  vector[N] R;  	
  vector[N] tem; 
  // priors
  vector[N] obs_error_i; 
  vector[N] nsamp; 
// priors
  real pr_C;
  real pr_C_sd;
  real pr_bpow;
  real pr_bpow_sd;
  real pr_bquad;
  real pr_bquad_sd;
  real pr_hyper_A;
  real pr_hyper_A_sd;
  real pr_hyper_A_deviation;
  real pr_hyper_A_deviation_sd;
  real pr_mean_A;
  real pr_mean_A_sd;
  real pr_hyper_B;
  real pr_hyper_B_sd;
  real pr_hyper_B_deviation;
  real pr_hyper_B_deviation_sd;
  real pr_mean_B;
  real pr_mean_B_sd;
  real pr_sigma_obs;
  real pr_sigma_obs_sd;
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
  real<lower=0> A;  
  real<lower=0> B; 
  real log_C;    
  real sigma_obs;
}
transformed parameters{
  vector[N] obs_error = (sigma_obs + obs_error_i)./nsamp; 
	real C=exp(log_C);
}
model {
  
  log_C ~ normal(pr_C, pr_C_sd);
  A ~ normal(pr_mean_A, pr_mean_A_sd);
  B ~ normal(pr_mean_B, pr_mean_B_sd);
  sigma_obs ~ normal(pr_sigma_obs,pr_sigma_obs_sd);
  
	for(i in 1:N)
	{
  y[i] ~ normal( A + B * exp(-t[i] * C), obs_error[i]);
    }
}
 generated quantities {
   vector[N] y_rep;
   vector[N] log_lik;
   
   for (i in 1:N){
     log_lik[i] = normal_lpdf(y[i] | A + B * exp(-t[i] * C), obs_error[i]);
  }   
   // Generate posterior predictive checks
   for (i in 1:N) {
	y_rep[i] = normal_rng(A + B * exp(-t[i] * C), obs_error[i]);
   } 
 }
 
