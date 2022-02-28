library(tidyverse)
library(R2jags)
library(bayesplot)
library(MCMCvis)
library(tidybayes)
library(ggpubr)
library(ggsci)
library(patchwork)

# Reading in the data

iso_data <- read.csv("../1 - Data/1 - Raw/Isotope Database Complete_14AUG20.csv")
dep_data <- read.csv("../1 - Data/1 - Raw/Departure Database Complete_14AUG20.csv")
migration_data <- readRDS("../1 - Data/1 - Raw/migration_rate.rds")
rainfall <- readRDS("../1 - Data/1 - Raw/rainfall.rds")
geos <- read.csv("../1 - Data/1 - Raw/AMRE_GEOS.csv")


rainfall_summary <- rainfall %>% 
  filter(B %in% c(1:3)) %>%
  group_by(YEAR) %>% summarize(total_rain = sum(RAIN..mm.),mean_rain = mean(RAIN..mm., na.rm=T), sd_rain = sd(RAIN..mm., na.rm=T) ) %>%
  mutate(cv=(sd_rain/mean_rain))  %>% data.frame()


dep_data$DEPDAT <- as.numeric(as.character(dep_data$DEPDAT))
dep_data$BAND <- as.factor(as.character(dep_data$BAND))
dep_data$HABITAT <- as.factor(dep_data$HABITAT)
dep_data$AGE <- as.factor(dep_data$AGE)
dep_data$SEX <- as.factor(dep_data$SEX)
levels(dep_data$SEX) <- c("Unknown","Female","Female","Male","Male","Female","Male")

levels(dep_data$AGE)[c(1,2,7,8,9,14,15)] <- c("Unknown")
levels(dep_data$AGE)[c(2:5)] <- c("Adult")
levels(dep_data$AGE)[c(3:6)] <- c("Juvenile")

migration_data <- migration_data %>% filter(SPECIES=="AMRE") %>% 
  select(DEPDAT2, mean_d2h, rate, SEX, HABITAT, mean_tailwind) %>% filter(!is.na(mean_d2h))

## Adding ID variable that includes year and band #

dep_data$ID <- paste(dep_data$BAND, dep_data$YEAR, sep="_")

iso_data$ID <- paste(iso_data$BAND, iso_data$YEAR, sep="_")

# Merging datasets and subsetting longitudinal data

compiled_data <- left_join(dep_data, iso_data, by="ID")

compiled_data <- compiled_data %>% select(YEAR.x,HABITAT, BAND.x, 
                                          AGE.x,SEX.x,DEPDAT,CENSOR,d2H,
                                          d2H.DUP,DUPVAR)
names(compiled_data) <- c("YEAR","HABITAT", "BAND","AGE","SEX",
                          "DEPDAT", "CENSOR","D2H","D2HDUP","D2HVAR")

repeats <- compiled_data %>% filter(AGE=="Adult",DEPDAT >= 20) %>% group_by(BAND) %>% mutate(RECAPS = n()) %>% arrange(desc(RECAPS))

repeats <- left_join(repeats,rainfall_summary, by = "YEAR") %>% filter(RECAPS > 1 & RECAPS < 19)

repeats <- repeats %>% group_by(BAND) %>% mutate(mD2H = mean(D2H, na.rm=T)) %>%
  filter(!is.nan(mD2H))

repeats$SEX <- as.factor(as.character(repeats$SEX))

ggplot(repeats, aes(x=scale(total_rain), y=DEPDAT, color=SEX)) +
  geom_point() + geom_smooth(method="lm", se=T) + scale_color_jco() +
  theme_pubr(base_size=20) + ylab("Departure Date \n(Days Since April 1st)") + 
  xlab("Scaled Winter Rainfall (Jan-Mar)") + theme(legend.title=element_blank())+
  facet_wrap(~HABITAT)


## Predicting Departure Date - Simple Random Intercepts Model

Nind <- length(levels(as.factor(repeats$BAND)))
levels(repeats$BAND) <- 1:Nind
repeats$total_rain_scaled <- scale(repeats$total_rain)[,1]

jagsdata <- with(repeats, list(departure = DEPDAT, rain = total_rain_scaled, D2H = mD2H, SEX = SEX,
                               indv = BAND, N = length(DEPDAT),HABITAT=HABITAT, Nind = Nind))

migration_data$SEX <- as.factor(migration_data$SEX)
levels(migration_data$SEX) <- c("Female","Male")
jags_migration_data <- with(migration_data[-c(6,12),], list(tailwind = mean_tailwind, d2H_m2 = mean_d2h, departure_m2 = DEPDAT2, rate = rate, SEXm = SEX, K = nrow(migration_data)-2))

jags_complete_data <- c(jagsdata,jags_migration_data)

jags_model <- function(){
  # Likelihood
  for (i in 1:N){
    departure[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + a[indv[i]] + beta * rain[i] + beta2 * SEX[i]  + beta4 * HABITAT[i] + beta3 * D2H[i] + beta5 * HABITAT[i]*rain[i]
  }
  
  # Priors for Model
  alpha ~ dnorm(0,100)
  beta ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
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
    fu[j] <- alphadelay + betadelay * delay[j] + betawind*tailwind[j] 
  }
  
  alphadelay ~ dnorm(0,100)
  betadelay ~ dnorm(0,100)
  betawind ~ dnorm(0,100)
  sigma2 ~ dunif(0,100)
  tau2 <- 1 / (sigma2*sigma2)
}


init_values <- function(){
  list(alpha = rnorm(1), sigma_a = runif(1), sigma2 = runif(1),
       beta = rnorm(1),beta2 = rnorm(1),beta3 = rnorm(1),beta4 = rnorm(1),beta5 = rnorm(1), sigma = runif(1),alphadelay = rnorm(1),betadelay = rnorm(1), betawind = rnorm(1))
}

params <- c("alpha","alphadelay", "beta", "beta2","beta3","beta4", "beta5","betadelay", "betawind","sigma", "sigma_a","vpc", "sigma2","delay")

fit_lm1 <- jags(data = jags_complete_data, inits = init_values, parameters.to.save = params, model.file = jags_model,
                n.chains = 3, n.iter = 100000, n.burnin = 5000, n.thin = 10, DIC = F)


MCMCsummary(fit_lm1, round = 2)

lm1_chains <- MCMCchains(fit_lm1, 
                         params = c("alpha","beta","beta2","beta3","beta4", "beta5"))