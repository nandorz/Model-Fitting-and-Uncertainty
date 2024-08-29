library(tidyverse)
library(gridExtra)

library(rstan)
library(bayesplot)
library(bridgesampling)


rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())


cases<- read.csv("case_data_set.csv")[,1]

N <- 20000; 

# times
n_days <- length(cases) 
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
cum0<- 1
M0 <- 0
y0 = c(S = s0, I = i0, R = r0, cumulativeI = cum0, M = M0)

# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases)

# number of MCMC steps
niter <- 3000


model1samples<-stan("model.stan",
                    iter=niter,
                    data=data_sir,
                    chains=3,
                    seed=123, 
                    control = list(max_treedepth = 13, adapt_delta = 0.99),
                    refresh=100,
                    cores = 3)


summary(model1samples)

model1samples

model_pred_SIR <- model1samples

pars=c('beta', 'gamma', "R0", "recovery_time", "alpha", "a","k")
print(model_pred_SIR, pars = pars)

#Checking stan outputs
beta_pred<-rstan::extract(model_pred_SIR)$beta
nu_pred<-rstan::extract(model_pred_SIR)$gamma
R0_pred<-rstan::extract(model_pred_SIR)$R0
inv_beta_pred<-rstan::extract(model_pred_SIR)$recovery_time
alpha_pred<-rstan::extract(model_pred_SIR)$alpha
a_pred<-rstan::extract(model_pred_SIR)$a
k_pred<-rstan::extract(model_pred_SIR)$k


# hist(beta_pred)
p1 <- mcmc_areas(model_pred_SIR,pars = "beta",prob = 0.95)    #95% credible interval 
p2 <- mcmc_areas(model_pred_SIR,pars = "gamma",prob = 0.95)     
p3 <- mcmc_areas(model_pred_SIR,pars = "R0",prob = 0.95)     
p4 <- mcmc_areas(model_pred_SIR,pars = "recovery_time",prob = 0.95)    
p5 <- mcmc_areas(model_pred_SIR,pars = "alpha",prob = 0.95) 
p6 <- mcmc_areas(model_pred_SIR,pars = "a",prob = 0.95) 
p7 <- mcmc_areas(model_pred_SIR,pars = "k",prob = 0.95) 
grid.arrange(p1,p2,p3,p4, p5, p6, p7)


#diagnostic plots
traceplot(model_pred_SIR, pars = c( "beta","gamma", "R0", "recovery_time","alpha", "a", "k"))
stan_dens(model_pred_SIR, pars =  c( "beta","gamma", "R0", "recovery_time", "alpha","a", "k"), separate_chains = TRUE)
pairs(model_pred_SIR, pars = c("beta", "gamma", "alpha", "a", "k"))

## Plot fitting results ####
# predicted cases ####
pred_cases_SIR<-rstan::extract(model_pred_SIR)$pred_cases
pred_cases_median<-apply(pred_cases_SIR, 2, median)
pred_cases_low<-apply(pred_cases_SIR,2, function(x) quantile(x,probs=0.025))
pred_cases_upp<-apply(pred_cases_SIR,2,function(x) quantile(x,probs=0.975))
# plot(pred_cases_upp,type="l",col="red")
# lines(pred_cases_median)
# lines(pred_cases_low,col="red")
# points(cases)

# ggplot
t_start <- 0
t_stop <- n_days
times <- seq(t_start, t_stop, by = 1)[-1]

data.Inc <-data.frame(days=times,
                      pred_median=pred_cases_median,
                      pred_low=pred_cases_low,
                      pred_upp=pred_cases_upp,
                      obs_cases=cases)

ggplot(data.Inc, aes(x=days, y = pred_median)) +
  geom_line(size = 0.5,color= "#00A5CF") +   
  geom_ribbon(aes(ymin=pred_low, ymax=pred_upp
  ), alpha=0.2, colour = NA,fill="#00A5CF")+
  geom_point( aes(days, obs_cases),color="#574AE2",size=1)+
  ylab("incidence") +
  xlab(" time")+
  ggtitle("Predicted Cases: The Information Model")

# Reff ####
Reff_SIR<-rstan::extract(model_pred_SIR)$Reff
Reff_median<-apply(Reff_SIR, 2, median)
Reff_low<-apply(Reff_SIR,2, function(x) quantile(x,probs=0.025))
Reff_upp<-apply(Reff_SIR,2,function(x) quantile(x,probs=0.975))
plot(Reff_upp,type="l",col="red")
lines(Reff_median)
lines(Reff_low,col="red")
#points(cases)

# ggplot
data<-data.frame(week=times,
                 pred_median=apply(Reff_SIR, 2, median),
                 pred_low=apply(Reff_SIR,2, function(x) quantile(x,probs=0.025)),
                 pred_upp=apply(Reff_SIR,2,function(x) quantile(x,probs=0.975)),
                 obs_cases=cases)

R0.modified.plot <- ggplot(data, aes(x=week, y = pred_median)) +
  geom_line(size = 0.5,color= "orange") +   
  geom_ribbon(aes(ymin=pred_low, ymax=pred_upp
  ), alpha=0.2, colour = NA,fill="orange")+
  #  geom_point( aes(week, obs_cases),color="#574AE2",size=1)+
  ylab("R0") +
  xlab(" time")+
  ggtitle("R0: The Information SIR Model")+
  theme_minimal()

R0.modified.plot

check_hmc_diagnostics(model_pred_SIR)

