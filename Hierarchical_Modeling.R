"
-----------------------------------------------------------
Fitting Hierarchical Bayesian Models to WM precision task


Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"
# Requirements before installing rstan (installed with brms) 
# https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows


if(!require('brms')) install.packages('brms')
if(!require('ggmcmc')) install.packages('ggmcmc')
if(!require('ggridges')) install.packages('ggridges')
if(!require('tidybayes')) install.packages('tidybayes')
if(!require('circular')) install.packages('circular')
if(!require('CircStats')) install.packages('CircStats')

library(brms)
library(ggridges)
library(tidyverse)
library(ggmcmc)
library(posterior)
library(bayesplot)
library(bayestestR)
library(tidybayes)

# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Data_Compiled/')
modelDir <- file.path(dataDir,'Models/')
plotDir <- file.path(parentDir,'Figures/')
setwd(dataDir)

# ----- Helper functions ------


# --- Read in Data ---

allSj_error.summary <- readRDS('loc_ss_respOrder_summary.RDS')

loc.ss_error.summary  <- readRDS('loc_ss_summary.RDS')

# have to change circular variables to numeric for WAIC comparisons
allSj_error.summary$circ_mean <- as.numeric(allSj_error.summary$circ_mean)
loc.ss_error.summary$circ_meanRad <- as.numeric(loc.ss_error.summary$circ_meanRad)

# ------------- Modeling ------------------

betaWeight_prior <- c(prior(normal(0,1)), class='b')


## MEAN OFFSET ##

"
Note, anovaBF does not test random factors; they are assumed to be nuisance factors.
Instead, it treats all the factors not specified as random as being fixed.
"

# **Need to specify that the prior for the betas come from a distribution that sums
# to zero, since ANOVA assumes that group deflections in the overall mean sum to 0.**
# get_prior(circ_meanRad ~ (1|Sj_Num), data=loc.ss_error.summary)

# Use von mises family for regression on circular data. Coefficients are expressed in radians

# intercept only model
mean.m0 <- brm(circ_meanRad ~ (1|Sj_Num), data=loc.ss_error.summary,
          family = von_mises,
          pior = betaWeight_prior,
          save_pars = save_pars(all=T), seed = 123,
          cores = 4, iter = 10000, chains = 4, warmup=1000,
          file=paste(modelDir,'meanOffset_intercept',sep='/'))

# The formula Response_Error ~ 1 + Set_Size + Location + (1 + Set_Size + Location|Sj_Num)
# Indicates that the fixed effects are set size and location, and that Individual Sjs vary randomly in terms of their intercept
# and their effect of Set Size and Location

mean.m1 <- brm(circ_meanRad ~ (1|Sj_Num) + SetSize, data=loc.ss_error.summary,
          family = von_mises,
          prior = betaWeight_prior,
          save_pars = save_pars(all=T), seed = 123,
          cores = 2, iter = 5000, chains = 4, warmup=1000,
          control = list(adapt_delta = 0.999,
                         max_treedepth = 13),
          file=paste(modelDir,'meanOffset_SS',sep='/'),
          file_refit = 'on_change')

mean.m2 <- brm(circ_meanRad ~ (1|Sj_Num) + Location, data=loc.ss_error.summary,
               family = von_mises,
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               control = list(adapt_delta = 0.9,
                              max_treedepth = 13),
               file=paste(modelDir,'meanOffset_Loc',sep='/'),
               file_refit = 'on_change')

mean.m3 <- brm(circ_meanRad ~ (1|Sj_Num) + SetSize + Location, data=loc.ss_error.summary,
          family = von_mises,
          prior = betaWeight_prior,
          save_pars = save_pars(all=T), seed = 123,
          cores = 2, iter = 5000, chains = 4, warmup=1000,
          control = list(adapt_delta = 0.999,
                         max_treedepth = 13),
          file=paste(modelDir,'meanOffset_Loc_SS',sep='/'),
          file_refit = 'on_change')

mean.m4 <- brm(circ_meanRad ~ (1|Sj_Num) + SetSize*Location, data=loc.ss_error.summary,
               family = von_mises,
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               control = list(adapt_delta = 0.999,
                              max_treedepth = 13),
               file=paste(modelDir,'meanOffset_LocXSS',sep='/'),
               file_refit = 'on_change')

# fit diagnostics
# Want to see fuzzy catapillars for the chains
plot(mean.m0)
plot(mean.m1)
plot(mean.m2)
plot(mean.m3)
plot(mean.m4)

# compare models using bayes factor

meanComp.bf <- bf_models(mean.m0,mean.m1,mean.m2,mean.m3,mean.m4,denominator = 1)
meanComp.loo <- LOO(mean.m0,mean.m1,mean.m2,mean.m3,mean.m4)


# ----- Summary of posterior distributions -------
post <- posterior_samples(m4, add_chain = T)


post %>% 
  ggplot(aes(x = b_Intercept)) +
  stat_halfeye(point_interval = mode_hdi, .width = c(.95)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(mu)) +
  theme_minimal_hgrid() +
  theme(plot.title.position = "plot")


# ------ Diagnostics --------
# MCMC should show good mixing and autocorrelation should die out after lag 1





## MRVL ##
"
Note, the range of MRVL can only be between 0 and 1, so might be best to use a beta regression
"
mrvl.m0 <- brm(MRVL ~ (1|Sj_Num), data = loc.ss_error.summary,
               family = 'beta',
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               file=paste(modelDir,'MRVL_intercept',sep='/'),
               control = list(adapt_delta = 0.999,
                              max_treedepth = 13))

mrvl.m1 <- brm(MRVL ~ SetSize + (1|Sj_Num), data = loc.ss_error.summary,
               family = 'beta',
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               file=paste(modelDir,'MRVL_SS',sep='/'),
               control = list(adapt_delta = 0.999,
                              max_treedepth = 13))

mrvl.m2 <- brm(MRVL ~ Location + (1|Sj_Num), data = loc.ss_error.summary,
               family = 'beta',
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               control = list(adapt_delta = 0.999,
                              max_treedepth = 13),
               file=paste(modelDir,'MRVL_Loc',sep='/'))

mrvl.m3 <- brm(MRVL ~ Location + SetSize + (1|Sj_Num), data = loc.ss_error.summary,
               family = 'beta',
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               control = list(adapt_delta = 0.999,
                              max_treedepth = 13),
               file=paste(modelDir,'MRVL_Loc_SS',sep='/'))

mrvl.m4 <- brm(MRVL ~ SetSize*Location + (1|Sj_Num), data = loc.ss_error.summary,
               family = 'beta',
               prior = betaWeight_prior,
               save_pars = save_pars(all=T), seed = 123,
               cores = 2, iter = 5000, chains = 4, warmup=1000,
               control = list(adapt_delta = 0.999,
                              max_treedepth = 13),
               file=paste(modelDir,'MRVL_LocXSS',sep='/'))

mrvlComp.bf <- bf_models(mrvl.m0,mrvl.m1,mrvl.m2,mrvl.m3,mrvl.m4,denominator = 1)
mrvlComp.loo <- LOO(mrvl.m0,mrvl.m1,mrvl.m2,mrvl.m3,mrvl.m4)

