functions {
  real switch_M(real alpha, real M) {
    return(1/(1 + alpha * M));
  }
  // defining the function and the input variables
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];           // defining that susceptible compartment (S) is the first value in vector y
      real I = y[2];           // defining that infected compartment (I) is the second value in vector y
      real R = y[3];           // defining that recovered compartment (R) is the third value in vector y
      real cumulativeI= y[4];  // defining a cumulative infected individuals (CumulativeI) is the fourth value in vector y
      real M = y[5];           // defining that index M is the fifth value in vector y
      
      real N = x_i[1];

      real beta = theta[1];    // beta, tansmission rate, is the first valye in vector theta
      real gamma = theta[2];   // gamma, recovery rate, is the second value in vector theta
      real alpha = theta[3];   // alpha, reactivity factor, is the third value in vector theta
      real a = theta[4];       // a, the inverse of average information delay, is the fourth value in vector theta
      real k = theta[5];       // k, the information coverage, is the fifth value in vector theta
      
      real dM_dt = a * (k * I - M); // this is how index M is defined

      real forcing_function = switch_M(alpha, M); // switch function for index M, it will be multiplied with beta at the next step
      
      real beta_eff = beta * forcing_function;    // incorporating information delay function into beta
      
      // the following are the equations for rates of change entering and leaving each compartment
      real dS_dt = -beta_eff * S / N;           
      real dI_dt =  beta_eff * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      real dcumulativeI_dt =  beta_eff * I * S / N ;
      return {dS_dt, dI_dt, dR_dt, dcumulativeI_dt, dM_dt};
  }
}
data {
  int<lower=1> n_days;
  real y0[5];
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
  // declaring that each of these parameters are real numbers whose value are equal to or more than 0

  real<lower=0> gamma;   // for gamma
  real<lower=0> beta;    // for beta
  real<lower=0> alpha;   // for alpha
  real<lower=0> a;       // for a
  real<lower=0> k;       // for k
  real<lower=0> phi_inv; // for phi_inv
}
transformed parameters{
  // declaring the value of each transformed parameter

  real<lower=0> y[n_days, 5]; // y is a matrix with 5 columns and n_day rows, containing real numbers more than or equal to 0
  real phi = 1. / phi_inv;    // defining phi as the inverse of phi_inv
  real new_infect[n_days];    // new_infect is a vector of containing real number as many as n_days 
  {
    // declaring that theta is a vector of three real numbers
    real theta[5];
    theta[1] = beta;  // assigning the first vector name as beta
    theta[2] = gamma; // assigning the second vector name as gamma 
    theta[3] = alpha; // assigning the third vector name as alpha
    theta[4] = a;     // assigning the fourth vector name as alpha
    theta[5] = k;     // assigning the fifth vector name as alpha

    y = integrate_ode_bdf(sir, y0, t0, ts, theta, x_r, x_i);
   // calculating the rate of change from S to I compartment
   new_infect[1] = to_matrix(y)[1,4]; // it extracts the first value of the cumulative incidence

   for (i in 2:n_days){
    new_infect[i] = to_matrix(y)[i,4]- to_matrix(y)[i-1,4]; // calculating the incidence by taking the difference between rows  
    }
  }
}

model {
  // setting up priors for each parameters
  beta ~ normal(2, 1);                   // prior for beta is a normal distribution with mean = 2, and sd = 1
  gamma ~ normal(0.1071429, 0.01785712); // as recovery rate (inverse of gamma) is between 7 and 14 days, so
                                         // the prior for gamma is a normal distibution with mean = (1/7 + 1/14)/2, and sd = (1/7 - 1/14)/2
  a ~ uniform(0,1);                      // prior for a has a uniform distribution from 0 to 1
  k ~ uniform(0,1);                      // prior for a has a uniform distribution from 0 to 1
  alpha ~ uniform(0,1);                  // prior for a has a uniform distribution from 0 to 1
  phi_inv ~ exponential(5);              // prior for a is an exponential distribution from 0 to 1

  for(j in 1:n_days){
     cases[j] ~ neg_binomial_2(0.1+new_infect[j], phi);
  }
}

generated quantities {
  real R0 = beta / gamma;
  real Reff[n_days];                // R0 but taking into account environmental changes
  real recovery_time = 1 / gamma;
  real inv_a = 1/a;
  real pred_cases[n_days];
  real M_t[n_days];
  real new_infections[n_days];
  new_infections[1] = to_matrix(y)[1,4];
  for (i in 2:n_days){
    new_infections[i] = to_matrix(y)[i,4]- to_matrix(y)[i-1,4];
  }
  for (i in 1:n_days){
    M_t[i] = y[i,5];
  }
  for (i in 1:n_days){
    Reff[i] = switch_M(alpha, M_t[i]) * beta / gamma;
  }
 pred_cases = neg_binomial_2_rng(new_infections, phi);
}


