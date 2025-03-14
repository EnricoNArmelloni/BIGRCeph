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
  real log_bF;   // Effect of covariate 1
  real log_bTB1;   // Effect of covariate 1
  real log_bTB2;   // Effect of covariate 1
  real b4TA;   // Effect of covariate 1
  vector [Ndiet-1] F_raw;
}
transformed parameters{
  vector[N] obs_error = (sigma_obs + obs_error_i)./nsamp; 
  real C=exp(log_C);
  real bF = exp(log_bF);
  real bTB1 = exp(log_bTB1);
  real bTB2 = -exp(log_bTB2);
  vector [Ndiet] F= inv_logit(append_row(F_raw, -sum(F_raw))); 
}
model {
  
  log_bF ~ normal(pr_bpow, pr_bpow_sd);
  log_bTB1 ~  normal(pr_bquad, pr_bquad_sd);
  log_bTB2 ~  normal(pr_bquad, pr_bquad_sd);
  bTA ~  normal(pr_bpow, pr_bpow_sd);
  log_C ~ normal(pr_C, pr_C_sd);
  
   A ~ normal(pr_mean_A, pr_mean_A_sd);
   B ~ normal(pr_mean_B, pr_mean_B_sd);
  
  F_raw ~ normal(pr_diet, pr_diet_sd);
  sigma_obs ~ normal(pr_sigma_obs,pr_sigma_obs_sd);
	for(i in 1:N)
	{
  y[i] ~ normal( A*((T_t[i]^bTA)+(R[i]*F[z[i]])^bF) + (B*(bTB1*T_t[i]+bTB2*T_t[i]^2 + (R[i]*F[z[i]])^bF)) * exp(-t[i] * C), obs_error[i]);
    }
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;

  for (i in 1:N){
    log_lik[i] = normal_lpdf(y[i] |  A*((T_t[i]^bTA)+(R[i]*F[z[i]])^bF) + (B*(bTB1*T_t[i]+bTB2*T_t[i]^2 + (R[i]*F[z[i]])^bF)) * exp(-t[i] * C), obs_error[i]);
 }
  
  // Generate posterior predictive checks
  for (i in 1:N) {
y_rep[i] = normal_rng( A*((T_t[i]^bTA)+(R[i]*F[z[i]])^bF) + (B*(bTB1*T_t[i]+bTB2*T_t[i]^2 + (R[i]*F[z[i]])^bF)) * exp(-t[i] * C), obs_error[i]);
  }
}
