"
-----------------------------------------------------------
Fitting Hierarchical Bayesian Model to WM precision task


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
if(!require('bpnreg')) install.packages('bpnreg')
library(brms)
library(tidyverse)
library(ggmcmc)
library(ggridges)
library(posterior)
library(bayesplot)
library(BayesFactor)
library(bayestestR)
library(tidybayes)
library(cowplot)
library(circular)
library(gridExtra)
library(grid)
library(ggpubr)
library(gtable)
library(bpnreg)

# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Data_Compiled/')
plotDir <- file.path(parentDir,'Figures/')
setwd(dataDir)

# ----- Helper functions ------

gamma_a_b_from_omega_sigma <- function(mode, sd) {
  if (mode <= 0) stop("mode must be > 0")
  if (sd   <= 0) stop("sd must be > 0")
  rate <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2)
  shape <- 1 + mode * rate
  return(list(shape = shape, rate = rate))
}


# --- Read in Data ---

# Mean Data
overall_loc_bias.df <- read.csv('AllSj_Overall_LocError_dataLong.csv', row.names = 1)
all_loc.by.ss_bias.df <- read.csv('AllSj_LocError_x_SS_dataLong.csv', row.names=1)

overall_loc_bias.df[c('Sj_num','Location')] <- lapply(overall_loc_bias.df[c('Sj_num','Location')], factor)
all_loc.by.ss_bias.df[c('Sj_num','Set.Size','Location')] <- lapply(all_loc.by.ss_bias.df[c('Sj_num','Set.Size','Location')],
                                                                   factor)


# Trial Data
long.allSj_locError <- read.csv('AllSj_Overall_LocError_Dist_dataLong.csv', row.names=1)
long.allSj_locError.ss_df <- read.csv('AllSj_Overall_LocError_x_Sss_Dist_dataLong.csv', row.names=1)

long.allSj_locError[c('Sj_Num','Location')] <- lapply(long.allSj_locError[c('Sj_Num','Location')], factor)
long.allSj_locError.ss_df[c('Sj_Num','Set_Size','Location')] <- lapply(long.allSj_locError.ss_df[c('Sj_Num','Set_Size','Location')],
                                                                   factor)

loc_bins <- seq(0,300,by=60)

all_sj.error_data <- readRDS('All_Sj_errorDist.RDS')

all_sj.raw_error <- all_sj.error_data$loc_overall$error

allSj_locError.raw_error.df <- lapply(all_sj.raw_error,function(x) as.data.frame(sapply(x,'[', 
                                                                                   seq(max(sapply(x, length))))))

allSj_locError.raw_error.df <- lapply(allSj_locError.raw_error.df, setNames,c('ResponseError')) 
                                        

for(iLoc in 1:length(allSj_locError.raw_error.df)){
  nTrials <- nrow(allSj_locError.raw_error.df[[iLoc]])
  loc_vec <- rep(loc_bins[iLoc],nTrials)
  allSj_locError.raw_error.df[[iLoc]]$Location <- loc_vec
  
  allSj_locError.raw_error.df[[iLoc]]$SampProbAng <- all_sj.error_data$loc_overall$samp_probeAngles[[iLoc]]
  allSj_locError.raw_error.df[[iLoc]]$RespProbAng <- all_sj.error_data$loc_overall$resp_probeAngles[[iLoc]]
  allSj_locError.raw_error.df[[iLoc]]$RespOrder <- all_sj.error_data$loc_overall$resp_order[[iLoc]]
  allSj_locError.raw_error.df[[iLoc]]$SetSize <- all_sj.error_data$loc_overall$setSizes[[iLoc]]
  allSj_locError.raw_error.df[[iLoc]]$Sj_Num <- all_sj.error_data$loc_overall$sj_nums[[iLoc]]
}

allSj_trialData.df <- bind_rows(allSj_locError.raw_error.df, .id='Location')

allSj_trialData.df <- allSj_trialData.df[c('Sj_Num','SampProbAng','RespProbAng', 'ResponseError',
                                           'Location','RespOrder','SetSize')]

allSj_trialData.df$RespOrder <- allSj_trialData.df$RespOrder + 1

allSj_trialData.df[c('Sj_Num','Location','SetSize','RespOrder')] <- lapply(allSj_trialData.df[c('Sj_Num',
                                                                                                'Location','SetSize',
                                                                                                'RespOrder')],
                                                                           factor)

levels(allSj_trialData.df$Location) <- loc_bins

# transform response errors greater than 180 for accurate angular difference
allSj_trialData.df$ResponseError <- unlist(lapply(allSj_trialData.df$ResponseError, function(x){
  
  if ( (x < 1) & (abs(x) > 180)){
    x <- (360- abs(x)) * -1 # keeps the direction of the angular difference. Clockwise = negative
  } else if ( (x > 1) & (abs(x) > 180)){
    x <- 360-abs(x)
  } else {
    x <- x
  }
  
}))

# for circular statistics functions
allSj_trialData.df$ResponseErrorRad <- allSj_trialData.df$ResponseError * pi/180 

rm(allSj_locError.raw_error.df,all_sj.error_data)

# K data
k_estimates <- read.csv('All_K_Estimates.csv')

subjects <- unique(allSj_trialData.df$Sj_Num)
# ------------- Aggregate Data -------------

allSj_error.summary <- allSj_trialData.df %>% 
  group_by(Sj_Num, Location, SetSize, RespOrder) %>% 
  summarise(
    circ_mean = round(mean_circ(ResponseErrorRad) * 180/pi,2),
    MRVL = round(rho_circ(ResponseErrorRad),2)
  )

loc.ss_error.summary <- allSj_trialData.df %>% 
  group_by(Sj_Num, Location, SetSize) %>% 
  summarise(
    circ_mean = round(mean_circ(ResponseErrorRad) * 180/pi,2),
    MRVL = round(rho_circ(ResponseErrorRad),2),
    sampProb_mean = round(mean_circ(SampProbAng * pi/180) * 180/pi,2),
    respProb_mean = round(mean_circ(RespProbAng * pi/180) * 180/pi,2)
  )



# -------------- EDA ----------------------

# Individual Response Error Distributions
ggdens_plot <- function (sj_num){
  data <- allSj_trialData.df[which(allSj_trialData.df$Sj_Num == sj_num),]
  
  ggplot(data, aes(x=ResponseError, color=Location))+
    geom_density(size=.5) +
    ylim(0,0.03) + 
    labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
         y = 'Density',
         title=paste('Sj: ',sj_num,'Overall')) +
    theme_classic() +
    theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
          axis.line=element_line(colour="black",size=.4),
          axis.text = element_text(size=7,color='black'),
          axis.title = element_text(size=8,color='black'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size=7))
}

p1 <- lapply(subjects,ggdens_plot)

ml <- marrangeGrob(p1,nrow=2,ncol=2)

ggsave("Individual_overallDists.pdf",plot=ml,path=plotDir)


# Individual ss plots
ggdens_plot.ss <- function (sj_num){
  data <- allSj_trialData.df[which(allSj_trialData.df$Sj_Num == sj_num),]
  
  ggplot(data, aes(x=ResponseError, color=Location))+
    geom_density(size=.5) + facet_wrap(SetSize ~ ., ncol=3) +
    ylim(0,0.06) + 
    labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
         y = 'Density',
         title=paste('Set Size Error','Sj: ',sj_num)) +
    theme_minimal() +
    theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
          axis.line=element_line(colour="black",size=.4),
          axis.text = element_text(size=7,color='black'),
          axis.title = element_text(size=8,color='black'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size=7))
}

p2 <- lapply(subjects,ggdens_plot.ss)

m2 <- marrangeGrob(p2,nrow=2,ncol=1)

ggsave("Individ_SS_Dists.pdf",plot=m2,path=plotDir)


# Mean Error
ggplot(data=loc.ss_error.summary, aes(x=Location, y=circ_mean, color=SetSize)) + 
  geom_point(size = 1.8) +
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("Response Offset ("*~degree*")",sep='')),
       title = 'Location by Set Size Mean Response Offset',
       fill = 'Set Size') +
  theme_classic() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

# Looking at frequency of response orders per location per set size
overall.resp_order <- ggplot(data=drop_na(allSj_error.summary), 
       aes(x=RespOrder, fill=Location, color=Location)) +
  geom_bar(stat='count',alpha=.5, position = position_dodge2()) +
  facet_wrap(~ SetSize) +
  labs(x = 'Response Order',
       y = 'Frequency',
       title = 'Frequency of Response Order\nFor Locations Across Set Sizes') +
  theme_minimal() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        strip.text = element_text(size=13, color = 'black'),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))
ggsave("Overall_RespOrder_SS_Dists.pdf",plot=overall.resp_order,path=plotDir)


# Individual response orders per location per set size
ggresp_hist <- function(sj_num){
  
  data <- allSj_trialData.df[which(allSj_trialData.df$Sj_Num == sj_num),]
  
  ggplot(data=drop_na(data), 
       aes(x=RespOrder, fill=Location, color=Location)) +
  geom_bar(stat='count',alpha=.5, position = position_dodge2()) +
  facet_wrap(~ SetSize, scales = 'free') +
  labs(x = 'Response Order',
       y = 'Frequency',
       title = paste('Sj: ', sj_num)) +
  theme_minimal() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        strip.text = element_text(size=13, color = 'black'),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))
  }

p3 <- lapply(subjects,ggresp_hist)

m3 <- marrangeGrob(p3,nrow=2,ncol=1)

ggsave("Individ_RespOrder_Freqs.pdf",plot=m3,path=plotDir)

# Raw Error as a function of set size, location, and response order
error.by.ssLocRespOrder <- ggplot(data=drop_na(allSj_trialData.df), 
                                  aes(x=ResponseError, color=Location)) + 
  geom_density(size = 1.1) + facet_grid(RespOrder ~ SetSize) +
  theme_minimal() + 
  labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
       y = 'Density',
       title='Response Error as a function of\nSet Size, Response Order, and Location') +
  theme_minimal() + 
  theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
      axis.line=element_line(colour="black",size=.4),
      axis.text = element_text(size=7,color='black'),
      axis.title = element_text(size=8,color='black'),
      legend.text = element_text(size = 6),
      legend.title = element_text(size=7))

ggsave('Overall_Error_by_SsLocRespOrder.jpg',plot=error.by.ssLocRespOrder,path=plotDir)

gg.respOrder.ss.loc <- function(sj_num){
  
  data <- allSj_trialData.df[which(allSj_trialData.df$Sj_Num == sj_num),]
  
  
  ggplot(data=drop_na(data), 
         aes(x=ResponseError, color=Location)) + 
    geom_density(size = 1.1) + facet_grid(RespOrder ~ SetSize) +
    labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
         y = 'Density',
         title=paste('Sj: ',sj_num)) +
    theme_minimal() + 
    theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
          axis.line=element_line(colour="black",size=.4),
          axis.text = element_text(size=7,color='black'),
          axis.title = element_text(size=8,color='black'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size=7))
}

p4 <- lapply(subjects,gg.respOrder.ss.loc)

m4 <- marrangeGrob(p4,nrow=1,ncol=1)

ggsave("Individ_RespOrder_ErrorDist.pdf",plot=m4,path=plotDir)

# ------------- Modeling ------------------

"
Note, anovaBF does not test random factors; they are assumed to be nuisance factors.
Instead, it treats all the factors not specified as random as being fixed.
"
# Replicate BF ANOVA
mean_offset <- mean(all_loc.by.ss_bias.df$Mean.Error)
sd_y <- sd(all_loc.by.ss_bias.df$Mean.Error)

omega <- sd_y / 2
sigma <- 2*sd_y

s_r <- gamma_a_b_from_omega_sigma(mode=omega, sd=sigma)


stanvars <- 
  stanvar(mean_offset, name = 'mean_offset') + 
  stanvar(sd_y, name = 'sd_y') + 
  stanvar(s_r$shape, name='alpha') + 
  stanvar(s_r$rate, name='beta')

# The formula Response_Error ~ 1 + Set_Size + Location + (1 + Set_Size + Location|Sj_Num)
# Indicates that the fixed effects are set size and location, and that Individual Sjs vary randomly in terms of their intercept
# and their effect of Set Size and Location
m0 <- brm(Mean.Error ~ 1,
          data=all_loc.by.ss_bias.df,
          family=gaussian,
          prior = c(prior(normal(mean_offset,sd_y*5), class=Intercept),
                    prior(gamma(alpha,beta), class=sd),
                    prior(cauchy(0, sd_y), class=sigma)),
          stanvars = stanvars,
          save_pars = save_pars(all=T),
          cores = 4)

m1 <- brm(Mean.Error ~ 1 + (1|Sj_num),
          data=all_loc.by.ss_bias.df,
          family=gaussian,
          prior = c(prior(normal(mean_offset,sd_y*5), class=Intercept),
                    prior(gamma(alpha,beta), class=sd),
                    prior(cauchy(0, sd_y), class=sigma)),
          stanvars = stanvars,
          save_pars = save_pars(all=T),
          cores = 4)

m2 <- brm(Mean.Error ~ 1 + (1 + Set.Size |Sj_num),
          data=all_loc.by.ss_bias.df,
          family=gaussian,
          prior = c(prior(normal(mean_offset,sd_y*5), class=Intercept),
                    prior(gamma(alpha,beta), class=sd),
                    prior(cauchy(0, sd_y), class=sigma)),
          stanvars = stanvars,
          save_pars = save_pars(all=T),
          cores = 4)



# compare models using bayes factor

bf <- bf_models(m0,m1,denominator = 1)

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

