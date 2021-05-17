"
-----------------------------------------------------------
Exploratory Data Analysis of WM precision task


Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"

library(tidyverse)
library(cowplot)
library(circular)
library(gridExtra)
library(grid)
library(ggpubr)
library(gtable)
library(CircStats)
library(rlist)
library(reshape2)

# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Data_Compiled/')
plotDir <- file.path(parentDir,'Figures/')
setwd(dataDir)


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
  allSj_locError.raw_error.df[[iLoc]]$DispCombo <- all_sj.error_data$loc_overall$combo_idx[[iLoc]]
}

allSj_trialData.df <- bind_rows(allSj_locError.raw_error.df, .id='Location')

allSj_trialData.df <- allSj_trialData.df[c('Sj_Num','SampProbAng','RespProbAng', 'ResponseError',
                                           'Location','RespOrder','SetSize','DispCombo')]

allSj_trialData.df$RespOrder <- allSj_trialData.df$RespOrder + 1

allSj_trialData.df[c('Sj_Num','Location','SetSize','RespOrder','DispCombo')] <- lapply(allSj_trialData.df[c('Sj_Num',
                                                                                                            'Location',
                                                                                                            'SetSize',
                                                                                                            'RespOrder',
                                                                                                            'DispCombo')],
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

allSj_trialData.df$SampProbRad <- circular(allSj_trialData.df$SampProbAng * pi/180,
                                           units ='rad')
allSj_trialData.df$RespProbRad <- circular(allSj_trialData.df$RespProbAng * pi/180,
                                           units='rad')

all_combos <- all_sj.error_data$loc_overall$combos
  
rm(allSj_locError.raw_error.df,all_sj.error_data)

# K data
k_estimates <- read.csv('All_K_Estimates.csv')

subjects <- unique(allSj_trialData.df$Sj_Num)

# ------------- Aggregate Data -------------

allSj_error.summary <- allSj_trialData.df %>% 
  group_by(Sj_Num, Location, SetSize, RespOrder) %>% 
  summarise(
    circ_mean = round(mean.circular(ResponseErrorRad) * 180/pi,2),
    MRVL = round(rho.circular(ResponseErrorRad),2)
  )

loc.ss_error.summary <- allSj_trialData.df %>% 
  group_by(Sj_Num, Location, SetSize) %>% 
  summarise(
    circ_mean = round(mean.circular(ResponseErrorRad) * 180/pi,2),
    MRVL = round(rho.circular(ResponseErrorRad),2),
    sampProb_mean = round(mean.circular(SampProbAng * pi/180) * 180/pi,2),
    respProb_mean = round(mean.circular(RespProbAng * pi/180) * 180/pi,2)
  )

# needed for modeling
loc.ss_error.summary$circ_meanRad <- loc.ss_error.summary$circ_mean * pi/180


saveRDS(allSj_error.summary,
        file = file.path(dataDir, 'loc_ss_respOrder_summary.RDS'))

saveRDS(loc.ss_error.summary,
        file = file.path(dataDir, 'loc_ss_summary.RDS'))

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


# -- Mean Direction x loc x SS --

meanOffset.loc_ss.plot <- ggplot(data=loc.ss_error.summary, 
       aes(x=Location, y=circ_mean, fill=SetSize)) +
  geom_hline(yintercept = 0, size = 1, colour='black') + 
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             shape=19, colour='slategrey', alpha = 0.6, show.legend=F) + 
  stat_summary(fun = function(x) mean.circular(x * pi/180) * 180/pi,
               geom='bar', position='dodge') +
  stat_summary(geom='errorbar', 
               fun.min = function(x) (mean.circular(x*pi/180)*180/pi) - (sd.circular(x*pi/180)*180/pi),
               fun.max = function(x) (mean.circular(x*pi/180)*180/pi) + (sd.circular(x*pi/180)*180/pi),
               position=position_dodge(width=1.1),width = .3, size=.8) +
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("Mean Response Offset ("*~degree*")",sep='')),
       title = 'Mean Response Offset by Location & Set Size',
       fill = 'Set Size') + 
  scale_x_discrete(labels = loc_bins) +
  scale_fill_manual(values = alpha(c('green4','deepskyblue3','red3'),.8)) + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'),
        panel.grid.major.x = element_line(size = 1, color='gray90'),
        axis.line.x = element_blank())

ggsave('loc_by_ss_meanError.jpg',
       plot = meanOffset.loc_ss.plot,
       path = plotDir)



# -- MRVL x loc x SS --
MRVL.loc_ss.plot <- ggplot(data=loc.ss_error.summary, 
                                 aes(x=Location, y=MRVL, fill=SetSize)) +
  stat_summary(fun = 'mean',
               geom='bar', position='dodge') +
  stat_summary(geom='errorbar', 
               fun.min = function(x) mean(x) - stats::sd(x)/sqrt(length(x)),
               fun.max = function(x) mean(x) + stats::sd(x)/sqrt(length(x)),
               position=position_dodge(width=1.1),width = .3, size=.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             shape=19, colour='black', show.legend=F) +
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = 'Resultant Vector Length',
       title = 'Resultant Vector Length by Location & Set Size',
       fill = 'Set Size') + 
  scale_x_discrete(labels = loc_bins) +
  scale_fill_manual(values = alpha(c('green4','deepskyblue3','red3'),.8)) + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('loc_by_ss_meanError.jpg',
       plot = MRVL.loc_ss.plot,
       path = plotDir)

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

# -------- Circular Correlation Matrices ----------------

allSj_loc.corrMatrix <- array(dim=c(length(subjects),
                                    2,
                                    length(loc_bins),
                                    length(loc_bins)))
for (iSub in 1:length(subjects)){
  sj_num <- subjects[iSub]
  
  for (iSs in 1:length(c(2,6))) {
    set_size <- c(2,6)[iSs]
    
    loc_combos <- all_combos[which(lapply(all_combos,
                                          function(x){length(x) == set_size}) == TRUE)]
    
    for (iComb in 1:length(loc_combos)){
      current_combo <- unlist(loc_combos[iComb])
      combo_idx <- list.findi(all_combos,identical(.,current_combo))
      
      for (iLoc in 1:length(current_combo)){
        
        samp_loc <- current_combo[iLoc]
        samp_loc.idx <- which(loc_bins == samp_loc)
        
        for (jLoc in 1:length(current_combo)) {
          resp_loc <- current_combo[jLoc]
          resp_loc.idx <- which(loc_bins == resp_loc)
          
          data <-  allSj_trialData.df[which(allSj_trialData.df$Sj_Num == sj_num & 
                                                allSj_trialData.df$DispCombo == combo_idx),]
          
          samp_probRad <- data$SampProbRad[which(data$Location == samp_loc)]
          resp_probRad <- data$RespProbRad[which(data$Location == resp_loc)]
          
          # handle missing values but keep trials matched
          if(any(is.na(resp_probRad))){
            na_idx <-  which(is.na(resp_probRad))
            samp_probRad <- samp_probRad[-na_idx]
            resp_probRad <- drop_na(resp_probRad)
          }
          
          rho <- cor.circular(samp_probRad,resp_probRad)
          
          allSj_loc.corrMatrix[iSub,iSs,samp_loc.idx,resp_loc.idx] <- round(rho,2)
          
        }
      }
    }
  }
}

# average over correlation matrices
avg.cor_mat <- apply(allSj_loc.corrMatrix, 2:4, mean)
avg.cor_mat.ss2 <- round(avg.cor_mat[1, ,],2)
avg.cor_mat.ss6 <- round(avg.cor_mat[2, ,],2)

rownames(avg.cor_mat.ss2) <- colnames(avg.cor_mat.ss2) <-  loc_bins
rownames(avg.cor_mat.ss6) <- colnames(avg.cor_mat.ss6) <-  loc_bins


# Plot matrices

melted_ss2 <- melt(avg.cor_mat.ss2)
melted_ss6 <- melt(avg.cor_mat.ss6)

melted_ss2$Var1 <- factor(melted_ss2$Var1)
melted_ss2$Var2 <- factor(melted_ss2$Var2)

melted_ss6$Var1 <- factor(melted_ss6$Var1)
melted_ss6$Var2 <- factor(melted_ss6$Var2)

# Var 2 should be the Sample Orientation Bin while Var 1 is the Response
colnames(melted_ss2) <- colnames(melted_ss6) <- c('Response',
                                                  'Sample',
                                                  'Rho')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                 "cyan", "#7FFF7F", "yellow", 
                                 "#FF7F00", "red", "#7F0000"))


ss2.cor_plot <- ggplot(melted_ss2, aes(x=Response, y=Sample, fill=Rho)) + 
  geom_tile(color = 'slategrey', size = .5)+
  geom_text(aes(label=Rho)) +
  scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-1,1)) +
  scale_y_discrete(limits = rev(levels(melted_ss2$Sample))) + 
  scale_x_discrete(position='top')+
  labs(x = expression(paste("Response Location ("*degree*")",sep='')), 
       y = expression(paste("Sample Location ("*~degree*")",sep='')),
       title='Set Size 2 Correlation Between \n Sample & Response Orientations\n Per Location') +
  theme_minimal() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  coord_fixed()

ggsave("ss2_corr_matrix.jpg",plot=ss2.cor_plot,path=plotDir)


ss6.cor_plot <- ggplot(melted_ss6, aes(x=Response, y=Sample, fill=Rho)) + 
  geom_tile(color = 'slategrey', size=.5)+
  geom_text(aes(label=Rho)) +
  scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-.5,.5)) +
  scale_y_discrete(limits = rev(levels(melted_ss2$Sample))) + 
  scale_x_discrete(position='top')+
  labs(x = expression(paste("Response Location ("*degree*")",sep='')), 
       y = expression(paste("Sample Location ("*~degree*")",sep='')),
       title='Set Size 6 Correlation Between \n Sample & Response Orientations\n Per Location') +
  theme_minimal() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  coord_fixed()
ggsave("ss6_corr_matrix.jpg",plot=ss6.cor_plot,path=plotDir)


