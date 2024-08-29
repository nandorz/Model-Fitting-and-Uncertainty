functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real cumulativeI= y[4]; 
      
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      real dcumulativeI_dt =  beta * I * S / N ;
      return {dS_dt, dI_dt, dR_dt, dcumulativeI_dt};
  }
}
data {
  int<lower=1> n_days;
  real y0[4];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}
transformed parameters{
  real<lower=0> y[n_days, 4];
  real phi = 1. / phi_inv;
  real new_infect[n_days];

  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
   new_infect[1] = to_matrix(y)[1,4];

   for (i in 2:n_days){
    new_infect[i] = to_matrix(y)[i,4]- to_matrix(y)[i-1,4];
  }
  }
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.1071429, 0.01785712);
  phi_inv ~ exponential(5);

  for(j in 1:n_days){
     cases[j] ~ neg_binomial_2(0.1+new_infect[j], phi);
  }
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  real new_infections[n_days];
  new_infections[1] = to_matrix(y)[1,4];
  for (i in 2:n_days){
    new_infections[i] = to_matrix(y)[i,4]- to_matrix(y)[i-1,4];
  }
 pred_cases = neg_binomial_2_rng(new_infections, phi);
}


