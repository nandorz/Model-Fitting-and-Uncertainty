functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real cumulativeI= y[4]; 

      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = 1/theta[2];
      
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
//  int cases[n_days];
}

transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1]=N;
}

parameters {
  real<lower=0> gamma_inv;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_days,4];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma_inv;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}


model {
  beta ~ normal(2, 1);
  gamma_inv ~ normal(10.5, 1.785714);
  phi_inv ~ exponential(5);

}

generated quantities {
  real R0 = beta*gamma_inv;
  real recovery_time = gamma_inv;
  real gamma = 1/ gamma_inv;
  real pred_cases[n_days];
  
  //col(matrix x, int n) - The n-th column of matrix x
  pred_cases = neg_binomial_2_rng(col(to_matrix(y),2) + 1e-5, phi);
}
