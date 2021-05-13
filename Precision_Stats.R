"
-----------------------------------------------------------
Stats on error estimates for WM precision task


Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"
if (!require('bayestestR')){install.packages('bayestestR')}
library(BayesFactor)
library(tidyverse)
library(bayestestR)
library(ggbeeswarm)
library(reshape)
library(ggplot2)
library(plotly)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)
library(gtable)

# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Data_Compiled/')
plotDir <- file.path(parentDir,'Figures/')
setwd(dataDir)

# read in data
all_sj.precision_data <- readRDS('All_Sj_Precision_Data.RDS')

all_sj.error_dist <- readRDS('All_Sj_errorDist.RDS')

k_estimates <- read.csv('All_K_Estimates.csv')

subjects <- all_sj.precision_data$sj_numbs

loc_bins <- seq(0,300,by=60)
# --------------- Helper Functions ---------------

se <- function(x){
  return(sd(x)/length(x))
}


BFanova.to.ggDf <- function(x){
  
  bfs <- extractBF(x)[1]
  model_names <- row.names(bfs)
  
  BF_anova.df <- as.data.frame(bfs)
  
  BF_anova.df$Models <- model_names
  colnames(BF_anova.df) <- c('BF','Model')
  
  return(BF_anova.df)
}

# jittering for polar plot
custom_point.jitter <- function(errors,angles,overlap_threshold = 1,jitter=5){
  checked_points <- list()
  jittered_points <- list()
  for (iError in 1:length(errors)){
    if (!is.element(iError,checked_points)){
      current_error <- errors[iError]
      closeness <- abs(current_error - errors)
      
      prob_points <- which(closeness < overlap_threshold & closeness != 0)
      
      prob_points <- setdiff(prob_points,checked_points)
      
      prob_points <- setdiff(prob_points,jittered_points)
      for (iPoint in prob_points){
        
        if ((!iPoint %% 2 ) & (iPoint > 1)){
          jitter <- jitter*2
        }
        
        angles[iPoint] <- angles[iPoint] + jitter
        
        jitter <- jitter * -1
      }
      
      jittered_points <- append(jittered_points,iPoint)
    }
    checked_points <- append(checked_points,iError)
    
  }
  
  return(angles)
}

# ----------- DF Creation -------------------------

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
all_loc.by.ss_bias.df <- data.frame()
for (iSub in 1:length(subjects)){
  sj_num <- subjects[iSub]
  loc.by.ss_mat_mean <- all_sj.precision_data$location.by.ss$overall$mean_error[iSub]
  loc.by.ss_mat_MRVL <- all_sj.precision_data$location.by.ss$overall$MRVL[iSub]
  
  # convert matrices into data frames
  sj_mean.df <- as.data.frame(loc.by.ss_mat_mean)
  sj_MRVL.df <- as.data.frame(loc.by.ss_mat_MRVL)
  
  names(sj_mean.df) <- names(sj_MRVL.df) <- seq(0,300,by=60)
  
  long_sj_mean.df <- melt(sj_mean.df, id.vars = NULL)
  long_sj_mean.MRVL <- melt(sj_MRVL.df, id.vars = NULL)
  
  sj_loc.ss.bias <- cbind(long_sj_mean.df,long_sj_mean.MRVL[2])
  
  colnames(sj_loc.ss.bias) <- c('Location','Mean.Error','MRVL')
  
  sj_loc.ss.bias$Location <- factor(sj_loc.ss.bias$Location, levels = seq(0,300,by=60))
  sj_loc.ss.bias$Set.Size <- factor(rep(c(1,2,6), nrow(sj_loc.ss.bias)/3))
  sj_loc.ss.bias$Sj_num <- factor(rep(sj_num,nrow(sj_loc.ss.bias)))
  
  sj_loc.ss.bias <- sj_loc.ss.bias[c('Sj_num','Set.Size','Location','Mean.Error','MRVL')]
  
  all_loc.by.ss_bias.df <- rbind(all_loc.by.ss_bias.df,sj_loc.ss.bias)
}

overall_loc_bias.df[c('Sj_num','Location')] <- lapply(overall_loc_bias.df[c('Sj_num','Location')], factor)

all_loc.by.ss_bias.df[c('Sj_num','Set.Size','Location')] <- lapply(all_loc.by.ss_bias.df[c('Sj_num','Set.Size',
                                                                                           'Location')], factor)
                                                                   

# save data frames to be used in other scripts
write.csv(overall_loc_bias.df,paste(dataDir,'AllSj_Overall_LocError_dataLong.csv', sep='/'))
write.csv(all_loc.by.ss_bias.df,paste(dataDir,'AllSj_LocError_x_SS_dataLong.csv', sep='/'))

# ------------------- Overall Error Stats ------------------------------------

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

# --- Model Plot ---
loc_bias.mean.bf_df <- BFanova.to.ggDf(loc_bias.mean.bfAnova)

loc_bias.mean.modelPlot <- ggplot(data = loc_bias.mean.bf_df , aes(x = BF, y = Model))+
  geom_bar(stat='identity', width = 0.5, color='black',
           fill='blue', lwd=.7, alpha = 0.3) +
  geom_text(aes(label=round(BF,2)), hjust = 1.5, 
            color = 'black', fontface='bold', size = 5)+
  ggtitle("vs Mean Error ~ Subject") + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('Loc_MeanError_Anova.jpg',
       loc_bias.mean.modelPlot,
       path=plotDir)

# MRVL
loc_bias.MRVL.bfAnova <- anovaBF(MRVL ~ Location + Sj_num, 
                            data=overall_loc_bias.df,
                            whichRandom = "Sj_num")

summary(loc_bias.MRVL.bfAnova)

# --- Model Plot ---
loc_bias.MRVL.bf_df <- BFanova.to.ggDf(loc_bias.MRVL.bfAnova)

loc_bias.MRVL.modelPlot <- ggplot(data = loc_bias.MRVL.bf_df , aes(x = BF, y = Model))+
  geom_bar(stat='identity', width = 0.5, color='black',
           fill='green', lwd=.7, alpha = 0.5) +
  geom_text(aes(label=round(BF,2)), hjust = 1.5, 
            color = 'black', fontface='bold', size = 5)+
  ggtitle("vs MRVL ~ Subject") + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('Loc_MRVL_Anova.jpg',
       loc_bias.MRVL.modelPlot,
       path=plotDir)

# ---- Mean Plots ----
loc_bias.mean.plot <- ggplot(data=loc_bias.summary, aes(x=Location, y=mean_error)) +
  geom_errorbar(aes(ymin = mean_error-mean_se, 
                    ymax = mean_error+mean_se), width = .1)+
  geom_point(size=5,fill='springgreen1',col='black',stroke=1,shape=22)+
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("Response Offset ("*~degree*")",sep='')),
       title = 'Location Mean Response Offset') +
  scale_x_discrete(labels = loc_bins) +
  geom_beeswarm(data=overall_loc_bias.df,aes(x=Location, y=Mean.Error),
             size=1.5,color='grey19',cex=3) +
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black')) 

ggsave('loc_overall_mean_error.jpg',
       plot = loc_bias.mean.plot,
       path = plotDir)

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

ggsave('loc_overall_MRVL.jpg',
       plot = loc_bias.MRVL.plot,
       path = plotDir)

# ------------------- Location x Set Size Stats ------------------------------------

# summary statistics
loc.ss.bias.summary <- all_loc.by.ss_bias.df %>% 
  group_by(Location, Set.Size) %>% 
  summarise(mean_error = mean(Mean.Error),
            mean_se = se(Mean.Error),
            mean_MRVL = mean(MRVL),
            MRVL_se = se(MRVL))


# Location Error for each Set Size
loc.ss.mean.bfAnova <- anovaBF(Mean.Error ~ Location + Set.Size + Sj_num, 
                                 data=all_loc.by.ss_bias.df,
                                 whichRandom = "Sj_num")
summary(loc.ss.mean.bfAnova)


loc.ss.mean.bf_df <- BFanova.to.ggDf(loc.ss.mean.bfAnova)

# ---- Model Plot -----
loc.ss.mean.modelPlot <- ggplot(data = loc.ss.mean.bf_df , aes(x = BF, y = Model))+
  geom_bar(stat='identity', width = 0.5, color='black',
           fill='blue', lwd=.7, alpha = 0.3) +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE)) +
  geom_text(aes(label=round(BF,2)), hjust = c(1.2,-0.1,-0.1,-0.1), 
            color = 'black', fontface='bold', size = 4)+
  ggtitle("vs Mean Error ~ Subject") + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('LocBySS_MeanError_Anova.jpg',
       loc.ss.mean.modelPlot,
       path=plotDir)


# Location MRVL for each Set Size
loc.ss.MRVL.bfAnova <- anovaBF(MRVL ~ Location + Set.Size + Sj_num, 
                               data=all_loc.by.ss_bias.df,
                               whichRandom = "Sj_num")

summary(loc.ss.MRVL.bfAnova)

loc.ss.MRVL.bf_df <- BFanova.to.ggDf(loc.ss.MRVL.bfAnova)

# --- Model Plot ---
loc.ss.MRVL.modelPlot <- ggplot(data = loc.ss.MRVL.bf_df , aes(x = BF, y = Model))+
  geom_bar(stat='identity', width = 0.5, color='black',
           fill='green', lwd=.7, alpha = 0.5) +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE)) +
  geom_text(aes(label=round(BF,2)), hjust = c(-0.1,-0.1,-0.1,1), 
            color = 'black', fontface='bold', size = 4)+
  ggtitle("vs MRVL ~ Subject") + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('LocBySS_MRVL_Anova.jpg',
       loc.ss.MRVL.modelPlot,
       path=plotDir)

# ---- Mean Plots ----
loc.ss_bias.mean.plot <- ggplot(data=loc.ss.bias.summary, 
                             aes(x=Location, y=mean_error, fill=Set.Size)) +
  geom_bar(position="dodge",stat="identity", colour = 'grey19') +
  geom_errorbar(aes(ymin = mean_error-mean_se, 
                    ymax = mean_error+mean_se), 
                position = position_dodge(width=0.9),width = .3, size=.8) +
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("Response Offset ("*~degree*")",sep='')),
       title = 'Location by Set Size Mean Response Offset',
       fill = 'Set Size') +
  scale_x_discrete(labels = loc_bins) +
  scale_fill_manual(values = alpha(c('green4','deepskyblue3','red3'),.7)) + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('loc_by_ss_meanError.jpg',
       plot = loc.ss_bias.mean.plot,
       path = plotDir)


loc.ss_bias.MRVL.plot <- ggplot(data=loc.ss.bias.summary, 
                             aes(x=Location, y=mean_MRVL, group=Set.Size)) +
  geom_point(aes(fill = Set.Size, shape = Set.Size), colour='black', size = 5, stroke=1.3) +
  geom_errorbar(aes(ymin = mean_MRVL-MRVL_se, 
                    ymax = mean_MRVL+MRVL_se),
                position= position_nudge(x=.15),
                width = .1, size = .8) +
  labs(x = expression(paste("Location ("*degree*")",sep='')), 
       y = expression(paste("MRVL",sep='')),
       title = 'Location by Set Size MRVL',
       shape = 'Set Size',
       fill = 'Set Size') +
  scale_x_discrete(labels = loc_bins) +
  scale_colour_manual(values = alpha(c('green4','deepskyblue3','red3'),.7)) + 
  scale_shape_manual(values = c(21,22,25)) + 
  theme_classic()+
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=1),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'))

ggsave('loc_by_ss_MRVL.jpg',
       plot = loc.ss_bias.MRVL.plot,
       path = plotDir)

# ------ Categorical Confusion Matrices ------
all.cc_matrices <- all_sj.precision_data$cc_matrix
all.cc_matrices <- array(unlist(all.cc_matrices),
                         c(length(subjects), 2, 6, 6))

# use dim 2:4 to preserve all but the subject dimension
avg.cc_matrices <- apply(all.cc_matrices, 2:4, mean)
ss2.cc_matrix <- avg.cc_matrices[1,,]
ss6.cc_matrix <- avg.cc_matrices[2,,]

rownames(ss2.cc_matrix) <- colnames(ss2.cc_matrix) <-  loc_bins
rownames(ss6.cc_matrix) <- colnames(ss6.cc_matrix) <-  loc_bins


# -- Plots --
melted_ss2 <- melt(ss2.cc_matrix)
melted_ss6 <- melt(ss6.cc_matrix)

melted_ss2$Var1 <- factor(melted_ss2$Var1)
melted_ss2$Var2 <- factor(melted_ss2$Var2)

melted_ss6$Var1 <- factor(melted_ss6$Var1)
melted_ss6$Var2 <- factor(melted_ss6$Var2)

# Var 2 should be the Sample Orientation Bin while Var 1 is the Response
colnames(melted_ss2) <- colnames(melted_ss6) <- c('Response',
                                                  'Sample',
                                                  'value')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                 "cyan", "#7FFF7F", "yellow", 
                                 "#FF7F00", "red", "#7F0000"))


ss2.cc_plot <- ggplot(melted_ss2, aes(x=Response, y=Sample, fill=value)) + 
  geom_tile(color = 'white')+
  scale_fill_gradientn(colors = jet.colors(10),
                       name='Frequency') +
  scale_y_discrete(limits = rev(levels(melted_ss2$Sample))) + 
  scale_x_discrete(position='top')+
  labs(x = expression(paste("Response Location ("*degree*")",sep='')), 
       y = expression(paste("Sample Location ("*~degree*")",sep='')),
       title='Set Size 2 Confusion Matrix') +
  theme_minimal() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  coord_fixed()

ggsave('ss2_confusion.jpg',
       plot = ss2.cc_plot,
       path = plotDir)


ss6.cc_plot <- ggplot(melted_ss6, aes(x=Response, y=Sample, fill=value)) + 
  geom_tile(color = 'white')+
  scale_fill_gradientn(colors = jet.colors(10),
                       name='Frequency') +
  scale_y_discrete(limits = rev(levels(melted_ss6$Sample))) + 
  scale_x_discrete(position='top')+
  labs(x = expression(paste("Response Location ("*degree*")",sep='')), 
       y = expression(paste("Sample Location ("*~degree*")",sep='')),
       title='Set Size 6 Confusion Matrix') +
  theme_minimal() +
  theme(plot.title = element_text(size=16,hjust = 0.5,face='bold'),
        axis.text = element_text(size=12,color='black'),
        axis.title = element_text(size=14,color='black'),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  coord_fixed()

ggsave('ss6_confusion.jpg',
       plot = ss6.cc_plot,
       path = plotDir)


# ------------ Response Error Distributions -------------
allSj_locError.overall <- all_sj.error_dist$loc_overall$error
allSj_locError.ss <- all_sj.error_dist$ss.loc$error

allSj_locError.overall <- data.frame(matrix(unlist(allSj_locError.overall),
                                            ncol = max(length(allSj_locError.overall)),
                                            byrow=TRUE))



colnames(allSj_locError.overall) <- loc_bins

long.allSj_locError <- melt(allSj_locError.overall)
colnames(long.allSj_locError) <- c('Location','Response_Error')
long.allSj_locError$Location <- factor(long.allSj_locError$Location,
                                       levels = loc_bins)

# include sj_nums to look at individual plots
long.allSj_locError$Sj_Num <- unlist(all_sj.error_dist$loc_overall$sj_nums, 
                                     recursive=FALSE)
long.allSj_locError$Sj_Num <- factor(long.allSj_locError$Sj_Num)

# lapply part will do it over each set size list of locations
# sapply creates a matrix and fills uneven rows with NA
allSj_locError.ss_dfs <- lapply(allSj_locError.ss,function(x) as.data.frame(sapply(x,'[', 
              seq(max(sapply(x, length))))))

# change column names
allSj_locError.ss_dfs <- lapply(allSj_locError.ss_dfs,setNames,loc_bins)

long.allSj_locError.ss_dfs <- lapply(allSj_locError.ss_dfs,melt)
long.allSj_locError.ss_dfs <- lapply(long.allSj_locError.ss_dfs,
                                     setNames,c('Location','Response_Error'))
long.allSj_locError.ss_df <- bind_rows(long.allSj_locError.ss_dfs,
                                       .id = 'Set_Size')
long.allSj_locError.ss_df$Set_Size <- as.factor(long.allSj_locError.ss_df$Set_Size)
levels(long.allSj_locError.ss_df$Set_Size) <- c(1,2,6)
long.allSj_locError.ss_df <- na.omit(long.allSj_locError.ss_df)


long.allSj_locError.ss_df$Sj_Num <- unlist(unlist(all_sj.error_dist$ss.loc$sj_nums, 
                                          recursive=F), 
                                   recursive=F)

long.allSj_locError.ss_df$Sj_Num <- factor(long.allSj_locError.ss_df$Sj_Num)

# save dataframes
write.csv(long.allSj_locError,paste(dataDir,'AllSj_Overall_LocError_Dist_dataLong.csv',sep='/'))
write.csv(long.allSj_locError.ss_df,paste(dataDir,'AllSj_Overall_LocError_x_Sss_Dist_dataLong.csv', sep='/'))


# --- Density Plots ---

# overall location error distribution
loc_overallDist.plot <- ggplot(long.allSj_locError, aes(x=Response_Error, color=Location))+
  geom_density(size=.7) +
  ylim(0,0.03) + 
  labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
       y = 'Density',
       title='Overall') +
  theme_minimal() +
  theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
        axis.line=element_line(colour="black",size=.4),
        axis.text = element_text(size=7,color='black'),
        axis.title = element_text(size=8,color='black'))


# location error x ss distributions
ss_dist_plts <- list()
for (iSet in 1:length(c(1,2,6))){
  
  current_ss <- c(1,2,6)[iSet]
  
  ss_dist_plts[[iSet]] <- local({
    
    iSet <- iSet
  
    p1 <- ggplot(long.allSj_locError.ss_df[which(long.allSj_locError.ss_df$Set_Size == current_ss),],
         aes(x=Response_Error, color=Location)) +
    geom_density(size = .7) +
      ylim(0,0.03) +
    labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
         y = 'Density',
         title=paste('Set Size:',current_ss, sep=' ')) +
    theme_minimal() +
    theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
          axis.line=element_line(colour="black",size=.4),
          axis.text = element_text(size=7,color='black'),
          axis.title = element_text(size=8,color='black'))
    
    if (iSet %in% c(1,2)){
      p1 <- p1 + theme(legend.position = 'none')
    } else {
      p1 <- p1 + theme(legend.text = element_text(size = 6),
                       legend.title = element_text(size=7)) + 
        scale_color_discrete(name = expression(paste("Location ("*degree*")", sep='')))
    }
    
    if (iSet %in% c(2,3)){
      p1 <- p1 + theme(
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
    }
    
    print(p1)
  })
}

plot_title <- textGrob('Location Response Error Distributions',
                       gp = gpar(fontsize = 15, fontface='bold'))

dist_grid <- grid.arrange(loc_overallDist.plot + theme(legend.position = 'none'), 
                          ss_dist_plts[[1]],
                          ss_dist_plts[[2]],
                          ss_dist_plts[[3]],
  top = plot_title,
  layout_matrix = rbind(c(NA,1,NA),
                        c(2,3,4))
  )

ggsave('Loc_Error_Dists.jpg',
       plot = dist_grid,
       path = plotDir)

## ---- Split as a function of K ----

high.K_sjs <- k_estimates$Sj_Numb[which(k_estimates$Kmax > median(k_estimates$Kmax))]
low.K_sjs <- k_estimates$Sj_Numb[which(k_estimates$Kmax <= median(k_estimates$Kmax))]

k.overall_plts <- list()
for(i in c(1:2)){
  if (i == 1){
    k_group <- low.K_sjs
    k_name <- 'Low K'
  } else {
    k_group <- high.K_sjs
    k_name <- 'High K'
  }
  
  k.overall_plts[[i]] <- local({
    
    i <- i
  
    plt <- ggplot(long.allSj_locError[which(long.allSj_locError$Sj_Num %in% k_group),],
                   aes(x=Response_Error, color=Location))+
      geom_density(size = .7) +
      ylim(0,0.025) + 
      labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
           y = 'Density',
           title=k_name) +
      theme_minimal() +
      theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
            axis.line=element_line(colour="black",size=.4),
            axis.text = element_text(size=7,color='black'),
            axis.title = element_text(size=8,color='black'))
  
    
    print(plt)
  })
  
}
plot_title <- textGrob('High and Low Capacity Subjects\nOverall Location Response Error',
                       gp = gpar(fontsize = 15, fontface='bold'))

gg_grid <- ggarrange(k.overall_plts[[1]],k.overall_plts[[2]], nrow = 1, 
                     common.legend = TRUE, legend='bottom')
k_loc_error_dists.plt  <- annotate_figure(gg_grid, top = plot_title)

ggsave('K_Loc_Error_Dists.jpg',
       plot = k_loc_error_dists.plt,
       path = plotDir)

# K x Set Size
k.ss_plots <- list()
plt_cnt <- 1
for (i in c(1:2)){
  
  if (i == 1){
    k_group <- low.K_sjs
    k_name <- 'Low K'
  } else {
    k_group <- high.K_sjs
    k_name <- 'High K'
  }

  for (iSet in 1:length(c(1,2,6))){
    
    current_ss <- c(1,2,6)[iSet]
    
    k.ss_plots[[plt_cnt]] <- local({
      
      plt_cnt <- plt_cnt 
      
      plt <- ggplot(long.allSj_locError.ss_df[which(long.allSj_locError.ss_df$Sj_Num %in% k_group & long.allSj_locError.ss_df$Set_Size == current_ss),],
                    aes(x=Response_Error, color=Location))+
        geom_density(size = .5) +
        labs(x = expression(paste("Response Error ("*degree*")",sep='')), 
             y = 'Density',
             title=paste('Set Size:',current_ss, sep=' ')) +
        ylim(c(0,0.04)) +
        theme_minimal() +
        theme(plot.title = element_text(size=9,hjust = 0.5,face='bold'),
              axis.line=element_line(colour="black",size=.4),
              axis.text = element_text(size=7,color='black'),
              axis.title = element_text(size=8,color='black'))
      if(plt_cnt < 6){
        plt <- plt + theme(legend.position = 'none')
      } else {
        plt <- plt + theme(legend.position = 'bottom', legend.title = element_text(size=6),
                           legend.text = element_text(size=5))
      }
    
      print(plt)
      
      
    })
    plt_cnt <- plt_cnt + 1
    
  }
}

plot_title <- textGrob('Location Set Size Response Error Distributions',
                       gp = gpar(fontsize = 15, fontface='bold'))

legend <- gtable_filter(ggplot_gtable(ggplot_build(k.ss_plots[[6]])), "guide-box")

low_title <- textGrob('Low K', gp = gpar(fontsize = 11, fontface='bold'))
high_title <- textGrob('High K', gp = gpar(fontsize = 11, fontface='bold'))

low_grid <- grid.arrange(k.ss_plots[[1]],k.ss_plots[[2]],k.ss_plots[[3]],
                         ncol = 3,
                         top = low_title)
high_grid <- grid.arrange(k.ss_plots[[4]],k.ss_plots[[5]],k.ss_plots[[6]] + theme(legend.position='none'),
                          ncol = 3,
                          top = high_title)

k_loc_ss_error_dists.plt <- grid.arrange(low_grid, high_grid, legend, top=plot_title, 
                                         heights=c(1, 1, 0.4))

ggsave('K_SS_Loc_Error_Dists.jpg',
       plot = k_loc_ss_error_dists.plt,
       path = plotDir)

# -- Individual subjects --

ggdens_plot <- function (sj_num){
  data <- long.allSj_locError[which(long.allSj_locError$Sj_Num == sj_num),]
  
  ggplot(data, aes(x=Response_Error, color=Location))+
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

# individual ss plots
ggdens_plot.ss <- function (sj_num){
  data <- long.allSj_locError.ss_df[which(long.allSj_locError.ss_df$Sj_Num == sj_num),]
  
  ggplot(data, aes(x=Response_Error, color=Location))+
    geom_density(size=.5) + facet_wrap(Set_Size ~ ., ncol=3) +
    ylim(0,0.03) + 
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

# ------ Polar Plots ----------------

# --- Overal Error ---

loc_bias.mean_error.polar <- plot_ly(
  data = loc_bias.summary,
  type='scatterpolar',
  mode = 'lines'
)

error <- loc_bias.summary %>% pull(mean_error)
error_se <- loc_bias.summary %>% pull(mean_se)

for(iLoc in 1:length(loc_bins)){
  offset <- error[iLoc]
  loc_angle <- loc_bins[iLoc] + 30
  values <- c(0,offset,offset,0)
  theta_angles <- c(0,loc_angle-5,loc_angle+5,0)
  
  # bar plot
  loc_bias.mean_error.polar <- loc_bias.mean_error.polar %>%
    add_trace(
      r = values,
      theta = theta_angles,
      fill = 'toself',
      fillcolor = 'rgba(0,0,255,0.5)',
      line = list(
        color = 'black'
      )
    )
  
  # sem
  offset_se <- error_se[iLoc]
  values <- c(offset-offset_se,offset+offset_se)
  theta_angles <- c(loc_angle,loc_angle)
  
  loc_bias.mean_error.polar <- loc_bias.mean_error.polar %>%
    add_trace(
      r=values,
      theta=theta_angles,
      line = list(
        color='black'
      )
    )
  
  # individual points
  all_offset <- overall_loc_bias.df$Mean.Error[which(overall_loc_bias.df$Location == iLoc)] 
  theta_angles <- rep(loc_angle,length(all_offset))
  jitter.theta_angles <- custom_point.jitter(all_offset,theta_angles,
                                             jitter = 4,
                                             overlap_threshold = 1)
  
  
  loc_bias.mean_error.polar <- loc_bias.mean_error.polar %>% 
    add_trace(
      r = all_offset,
      theta = jitter.theta_angles,
      line = list(
        width = 0
      ),
      marker = list(
        color = 'black',
        symbol = 'circle',
        size = 5,
        opacity = 0.8
      )
    )
}

loc_bias.mean_error.polar <- loc_bias.mean_error.polar %>% 
  layout(
  title = list(
    text = '<b>Location Mean Response Degree Offset</b>',
    font = list(
      size = 16,
      color='black'
    ),
    x = 0.54),
  polar = list(
    angularaxis = list(
      tickmode = 'array',
      tickvals = loc_bins,
      ticktext = sapply(loc_bins,
                        function(x){paste(as.character(x),'°',sep = '')}),
      tickwidth = 2,
      linewidth = 2,
      layer = 'below traces',
      tickfont = list(
        size = 14 
      )
    ),
    radialaxis = list(
      dtick=10,
      tick0 = 0,
      showline = F,
      tickfont = list(
        size = 12
      )
    )
  ),
  showlegend=F
)

loc_bias.mean_error.polar


# --- Overall MRVL ---

loc_bias.MRVL.polar <- plot_ly(
  data = loc_bias.summary,
  type='scatterpolar',
  mode = 'lines'
)

MRVL <- loc_bias.summary %>% pull(mean_MRVL)
MRVL_se <- loc_bias.summary %>% pull(MRVL_se)

for(iLoc in 1:length(loc_bins)){
  loc_mrvl <- MRVL[iLoc]
  loc_angle <- loc_bins[iLoc] + 30
  values <- c(0,loc_mrvl,loc_mrvl,0)
  theta_angles <- c(0,loc_angle-5,loc_angle+5,0)
  
  # bar plot
  loc_bias.MRVL.polar <- loc_bias.MRVL.polar %>%
    add_trace(
      r = values,
      theta = theta_angles,
      fill = 'toself',
      fillcolor = 'rgba(0,255,0,0.5)',
      line = list(
        color = 'black'
      )
    )
  
  # sem
  loc.MRVL_se <- MRVL_se[iLoc]
  values <- c(loc_mrvl-loc.MRVL_se,loc_mrvl+loc.MRVL_se)
  theta_angles <- c(loc_angle,loc_angle)
  
  loc_bias.MRVL.polar <- loc_bias.MRVL.polar %>%
    add_trace(
      r=values,
      theta=theta_angles,
      line = list(
        color='black'
      )
    )
  
  # individual points
  all_MRVL <- overall_loc_bias.df$MRVL[which(overall_loc_bias.df$Location == iLoc)] 
  theta_angles <- rep(loc_angle,length(all_MRVL))
  jitter.theta_angles <- custom_point.jitter(all_MRVL,theta_angles,
                                             jitter = 1,
                                             overlap_threshold = 0.1)
  
  
  loc_bias.MRVL.polar <- loc_bias.MRVL.polar %>% 
    add_trace(
      r = all_MRVL,
      theta = jitter.theta_angles,
      line = list(
        width = 0
      ),
      marker = list(
        color = 'black',
        symbol = 'circle',
        size = 5,
        opacity = 0.8
      )
    )
}

loc_bias.MRVL.polar <- loc_bias.MRVL.polar %>% 
  layout(
    title = list(
      text = '<b>Location MRVL</b>',
      font = list(
        size = 16,
        color='black'
      ),
      x = 0.54),
    polar = list(
      angularaxis = list(
        tickmode = 'array',
        tickvals = loc_bins,
        ticktext = sapply(loc_bins,
                          function(x){paste(as.character(x),'°',sep = '')}),
        tickwidth = 2,
        linewidth = 2,
        layer = 'below traces',
        tickfont = list(
          size = 14 
        )
      ),
      radialaxis = list(
        dtick=.2,
        tick0 = 0,
        showline = F,
        tickfont = list(
          size = 12
        )
      )
    ),
    showlegend=F
  )

loc_bias.MRVL.polar



