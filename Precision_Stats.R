"
-----------------------------------------------------------
Stats on error estimates for WM precision task


Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"
if (!require('bayestestR')) install.packages('bayestestR')
library(BayesFactor)
library(tidyverse)
library(bayestestR)
library(ggbeeswarm)


# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Data_Compiled/')

setwd(dataDir)

# read in data
all_sj.precision_data <- readRDS('All_Sj_Precision_Data.RDS')

subjects <- all_sj.precision_data$sj_numbs

loc_bins <- seq(0,300,by=60)
# --------------- Helper Functions ---------------

se <- function(x){
  return(sd(x)/length(x))
}


# ----------- Stats -------------------------

# Test for overall bias across locations
overall_loc_bias.df <- data.frame()
all_sj.loc_meanError <- all_sj.precision_data$location_overall$mean_error
all_sj.loc_MRVL <- all_sj.precision_data$location_overall$MRVL

row_count <- 1
for(iSub in 1:length(subjects)){
  sj_num <- subjects[iSub]
  sj.loc_mean_error <- unlist(all_sj.loc_meanError[iSub])
  sj.loc_MRVL <- unlist(all_sj.loc_MRVL[iSub])
  for (iLoc in c(1:6)){
    overall_loc_bias.df[row_count,1] <- sj_num
    overall_loc_bias.df[row_count,2] <- iLoc
    overall_loc_bias.df[row_count,3] <- sj.loc_mean_error[iLoc]
    overall_loc_bias.df[row_count,4] <- sj.loc_MRVL[iLoc]
    
    row_count <- row_count + 1
  }
  
}
colnames(overall_loc_bias.df) <- c('Sj_num','Location','Mean.Error','MRVL')
overall_loc_bias.df$Sj_num <- factor(overall_loc_bias.df$Sj_num, levels = subjects)
overall_loc_bias.df$Location <- factor(overall_loc_bias.df$Location, levels = c(1:6))


# Test for location bias as a function of set size
loc.by.ss_bias.df <- data.frame()
row_count <- 1
for(iSub in 1:length(subjects)){
  sj_num <- subjects[iSub]
  
  sj_ss1.loc.error <- unlist(all_sj.precision_data$location.by.ss$ss1[iSub])
  sj_ss2.loc.error <- unlist(all_sj.precision_data$location.by.ss$ss2[iSub])
  sj_ss6.loc.error <- unlist(all_sj.precision_data$location.by.ss$ss6[iSub])
  for (iLoc in c(1:6)){
     # CAN JUST USE DIAGONAL FOR SS 1?
    
  }
  
}

# summary statistics
loc_bias.summary <- overall_loc_bias.df %>% 
  group_by(Location) %>% 
  summarise(mean_error = mean(Mean.Error),
            mean_se = se(Mean.Error),
            mean_MRVL = mean(MRVL),
            MRVL_se = se(MRVL))

# Error
loc_bias.mean.bfAnova <- anovaBF(Mean.Error ~ Location + Sj_num, 
                            data=overall_loc_bias.df,
                            whichRandom = "Sj_num")
summary(loc_bias.mean.bfAnova)

par(ask=F)
plot(loc_bias.mean.bfAnova,col='light blue',ylab='BF')


loc_bias.mean.bfAnova_chains <- posterior(loc_bias.mean.bfAnova, iterations = 10000)
summary(loc_bias.mean.bfAnova_chains[,1:7])

plot(loc_bias.mean.bfAnova_chains[,2:7])

ci(loc_bias.mean.bfAnova_chains[,1:7])

# MRVL
loc_bias.MRVL.bfAnova <- anovaBF(MRVL ~ Location + Sj_num, 
                            data=overall_loc_bias.df,
                            whichRandom = "Sj_num")

summary(loc_bias.MRVL.bfAnova)

par(ask=F)
plot(loc_bias.MRVL.bfAnova,col='light green',ylab='BF')

loc_bias.MRVL.bfAnova_chains <- posterior(loc_bias.MRVL.bfAnova, iterations = 10000)
summary(loc_bias.MRVL.bfAnova_chains[,1:7])

plot(loc_bias.MRVL.bfAnova_chains[,2:7])

ci(loc_bias.MRVL.bfAnova_chains[,1:7])


# plot mean error and MRVL at each location
loc_bias.mean.plot <- ggplot(data=loc_bias.summary, aes(x=Location, y=mean_error)) +
  geom_errorbar(aes(ymin = mean_error-mean_se, 
                    ymax = mean_error+mean_se), width = .1)+
  geom_point(size=5,fill='springgreen1',col='black',stroke=1,shape=22)+
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("Response Offset ("*~degree*")",sep='')),
       title = 'Location Mean Response Offset') +
  scale_x_discrete(labels = loc_bins) +
  geom_beeswarm(data=overall_loc_bias.df,aes(x=Location, y=Mean.Error),
             size=1.5,color='slategrey',cex=3) +
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))
  

loc_bias.MRVL.plot <- ggplot(data=loc_bias.summary, aes(x=Location, y=mean_MRVL)) +
  geom_errorbar(aes(ymin = mean_MRVL-MRVL_se, 
                    ymax = mean_MRVL+MRVL_se), width = .1)+
  geom_point(size=5,fill='tomato3',col='black',stroke=1,shape=22)+
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("MRVL",sep='')),
       title = 'Location MRVL') +
  scale_x_discrete(labels = loc_bins) +
  geom_beeswarm(data=overall_loc_bias.df,aes(x=Location, y=MRVL),
                size=1.5,color='slategrey',cex=3) +
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))
  