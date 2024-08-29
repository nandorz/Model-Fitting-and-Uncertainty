library(tidyverse)
library(gridExtra)

library(rstan)
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
y0 = c(S = s0, I = i0, R = r0, cumulativeI = cum0)

# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N)#, cases = cases)

# number of MCMC steps
niter <- 3000


model1samples.prior1 <-stan("sir_prior.stan",
                    iter=niter,
                    data=data_sir,
                    chains=3,
                    seed=123, 
                    refresh=100)


summary(model1samples.prior1)


model1samples.prior1

model_pred_SIR.prior1 <- model1samples.prior1

pars=c('beta', 'gamma', "R0", "recovery_time")
print(model_pred_SIR.prior1, pars = pars)

#Checking stan outputs
beta_pred<-rstan::extract(model_pred_SIR.prior1)$beta
nu_pred<-rstan::extract(model_pred_SIR.prior1)$gamma
R0_pred<-rstan::extract(model_pred_SIR.prior1)$R0
inv_beta_pred<-rstan::extract(model_pred_SIR.prior1)$recovery_time


# hist(beta_pred)
p1 <- mcmc_areas(model_pred_SIR.prior1,pars = "beta",prob = 0.95)    #95% credible interval 
p2 <- mcmc_areas(model_pred_SIR.prior1,pars = "gamma",prob = 0.95)     
p3 <- mcmc_areas(model_pred_SIR.prior1,pars = "R0",prob = 0.95)     
p4 <- mcmc_areas(model_pred_SIR.prior1,pars = "recovery_time",prob = 0.95)    
grid.arrange(p1,p2,p3,p4)


#diagnostic plots
traceplot(model_pred_SIR.prior1, pars = c( "beta","gamma", "R0", "recovery_time"))
stan_dens(model_pred_SIR.prior1, pars =  c( "beta","gamma", "R0", "recovery_time"), separate_chains = TRUE)
pairs(model_pred_SIR.prior1, pars = c("beta", "gamma"))


#Plot fitting results
pred_cases_SIR<-rstan::extract(model_pred_SIR.prior1)$pred_cases
pred_cases_median<-apply(pred_cases_SIR, 2, median)
pred_cases_low<-apply(pred_cases_SIR,2, function(x) quantile(x,probs=0.025))
pred_cases_upp<-apply(pred_cases_SIR,2,function(x) quantile(x,probs=0.975))
plot(pred_cases_upp,type="l",col="red")
lines(pred_cases_median)
lines(pred_cases_low,col="red")
points(cases)

#ggplot
t_start <- 0
t_stop <- n_days
times <- seq(t_start, t_stop, by = 1)[-1]


data<-data.frame(week=times,
                 pred_median=apply(pred_cases_SIR, 2, median),
                 pred_low=apply(pred_cases_SIR,2, function(x) quantile(x,probs=0.025)),
                 pred_upp=apply(pred_cases_SIR,2,function(x) quantile(x,probs=0.975)),
                 obs_cases=cases)

ggplot(data, aes(x=week, y = pred_median)) +
  geom_line(size = 0.5,color= "#00A5CF") +   
  geom_ribbon(aes(ymin=pred_low, ymax=pred_upp
  ), alpha=0.2, colour = NA,fill="#00A5CF")+
  geom_point( aes(week, obs_cases),color="#574AE2",size=1)+
  ylab("incidence") +
  xlab(" time")+
  ggtitle("Model fitting based on SIR: The Information Model")

# Reff ####
Reff_SIR<-rstan::extract(model_pred_SIR.prior1)$Reff
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

ggplot(data, aes(x=week, y = pred_median)) +
  geom_line(size = 0.5,color= "#00A5CF") +   
  geom_ribbon(aes(ymin=pred_low, ymax=pred_upp
  ), alpha=0.2, colour = NA,fill="#00A5CF")+
  #  geom_point( aes(week, obs_cases),color="#574AE2",size=1)+
  ylab("R0") +
  xlab(" time")+
  ggtitle("R0: The Information Model")

check_hmc_diagnostics(model_pred_SIR.prior1)

# prior predictive check ####
s_prior <- rstan::extract(model_pred_SIR.prior1)
phi_inv_prior.check <- ggplot(tibble(r = s_prior$phi_inv)) +
  geom_density(aes(x = r), fill = "#00A5CF", alpha = 0.6) +
  geom_vline(xintercept = c(0.01 ,20), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Dispersion parameter (log)", y = "Probability density")

beta_prior.check <- ggplot(tibble(r = s_prior$beta)) +
  geom_density(aes(x = r), fill = "#00A5CF", alpha = 0.6) +
  geom_vline(xintercept = c(0.01 ,20), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Transmission rate (log)", y = "Probability density")

gamma_prior.check <- ggplot(tibble(r = s_prior$gamma)) +
  geom_density(aes(x = r), fill = "#00A5CF", alpha = 0.6) +
  geom_vline(xintercept = c(0.01 ,20), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Recovery rate (log)", y = "Probability density")

recovery_time_prior.check <- ggplot(tibble(r = s_prior$recovery_time)) +
  geom_density(aes(x = r), fill = "#00A5CF", alpha = 0.6) +
  geom_vline(xintercept = c(0.01 ,20), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Recovery time (days, log)", y = "Probability density")

R0_prior.check <- ggplot(tibble(r = s_prior$R0)) +
  geom_density(aes(x = r), fill = "#00A5CF", alpha = 0.6) +
  geom_vline(xintercept = c(0.01 ,20), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Basic reproduction number,R0 , (log)", y = "Probability density")

grid.arrange(beta_prior.check, gamma_prior.check, recovery_time_prior.check, R0_prior.check, phi_inv_prior.check)

n_draws <- 1000
draws <- as_tibble(t(s_prior$y[,,2][1:n_draws,])) %>% add_column(t=t)
draws <- pivot_longer(draws, c(1:1000), names_to = "draw")
draws %>%
  ggplot() +
  geom_line(mapping = aes(x = t, y = value, group = draw), alpha = 0.6, size = 0.1) +
  geom_hline(yintercept = 20000, color = "red") +
  geom_text(x = 100, y = 19000, label = "Population size", color = "red") +
  labs(x= "Day", y = "Number of infected students")


smr_pred <- cbind(as.data.frame(summary(model_pred_SIR.prior1, pars = "pred_cases",
                                        probs = c(0.05, 0.5, 0.95))$summary), t)

colnames(smr_pred) <- make.names(colnames(smr_pred))


ggplot(smr_pred, mapping = aes(x=t)) + 
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "#00A5CF", alpha = 0.35) +
  geom_line(mapping = aes(x=t, y = X50.), color = "#00A5CF") +
  geom_hline(yintercept = 20000, color = "red") +
  geom_text(x= 100, y = 19600, label = "Population size", color = "red") +
  labs(x = "Day", y = "Number of students in bed")