library(tidyverse)
library(R2jags)
library(bayesplot)
library(MCMCvis)
library(tidybayes)
library(ggpubr)
library(ggsci)
library(patchwork)

# Reading in the data

jags_complete_data <- readRDS("./jags_data.rds")

jags_model <- function(){
  # Likelihood
  for (i in 1:N){
    departure[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + a[indv[i]] + beta * rain[i] + beta2 * SEX[i]  + beta4 * HABITAT[i] + beta3 * D2H[i] + beta5 * rain[i] * HABITAT[i]
  }
  
  # Priors for Model
  alpha ~ dnorm(0,0.01)
  beta ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  beta3 ~ dnorm(0,0.01)
  beta4 ~ dnorm(0,0.01)
  beta5 ~ dnorm(0,0.01)
  sigma ~ dunif(0,100)
  tau <- 1 / (sigma*sigma)
  
  # Priors for Random Intercept
  sigma_a ~ dunif(0,100)
  tau_a <- 1/(sigma_a * sigma_a)
  for (j in 1:Nind){
    a[j] ~ dnorm(0, tau_a)
  }
  
  vpc <- sigma / (sigma + sigma_a)
  
  for (j in 1:K){
    rate[j] ~ dnorm(fu[j], tau2)
    delay[j] <- departure_m2[j] - (alpha+beta2*SEXm[j] + beta3*d2H_m2[j])
    fu[j] <- alphadelay + betadelay * delay[j] + betawind * tailwind[j]
  }
  
  alphadelay ~ dnorm(0,0.01)
  betadelay ~ dnorm(0,0.01)
  betawind ~ dnorm(0,0.01)
  sigma2 ~ dunif(0,100)
  tau2 <- 1 / (sigma2*sigma2)
}


init_values <- function(){
  list(alpha = rnorm(1), sigma_a = runif(1), sigma2 = runif(1),
       beta = rnorm(1),beta2 = rnorm(1),beta3 = rnorm(1),beta4 = rnorm(1), beta5 = rnorm(1), sigma = runif(1),alphadelay = rnorm(1),betadelay = rnorm(1),betawind = rnorm(1))
}

params <- c("alpha","alphadelay", "beta", "beta2","beta3","beta4","beta5","betadelay","betawind", "sigma", "sigma_a","vpc", "sigma2","delay")

fit_lm1 <- jags(data = jags_complete_data, inits = init_values, parameters.to.save = params, model.file = jags_model,
                n.chains = 6, n.iter = 100000, n.burnin = 5000, n.thin = 10, DIC = F)


MCMCsummary(fit_lm1, round = 3)
