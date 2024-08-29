library(tidyverse)
library(gridExtra)

library(rstan)

rm(list = ls())
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


model1samples.prior2 <-stan("sir_prior2.stan",
                    iter=niter,
                    data=data_sir,
                    chains=3,
                    seed=123, 
                    refresh=100)


summary(model1samples.prior2)


model1samples.prior2

model_pred_SIR.prior2 <- model1samples.prior2

pars=c('beta', 'gamma', "R0", "recovery_time")
print(model_pred_SIR.prior2, pars = pars)

#Checking stan outputs
beta_pred<-rstan::extract(model_pred_SIR.prior2)$beta
nu_pred<-rstan::extract(model_pred_SIR.prior2)$gamma
R0_pred<-rstan::extract(model_pred_SIR.prior2)$R0
inv_beta_pred<-rstan::extract(model_pred_SIR.prior2)$recovery_time

# prior predictive check ####
s_prior <- rstan::extract(model_pred_SIR.prior2)
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


smr_pred <- cbind(as.data.frame(summary(model_pred_SIR.prior2, pars = "pred_cases",
                                        probs = c(0.05, 0.5, 0.95))$summary), t)

colnames(smr_pred) <- make.names(colnames(smr_pred))


ggplot(smr_pred, mapping = aes(x=t)) + 
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "#00A5CF", alpha = 0.35) +
  geom_line(mapping = aes(x=t, y = X50.), color = "#00A5CF") +
  geom_hline(yintercept = 20000, color = "red") +
  geom_text(x= 100, y = 19600, label = "Population size", color = "red") +
  labs(x = "Day", y = "Number of students in bed")