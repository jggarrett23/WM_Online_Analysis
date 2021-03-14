"
-----------------------------------------------------------
Compute error estimates for WM precision task


Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"

# Load necessary packages
if (!require('rjson')) install.packages('rjson')
if (!require('googlesheets4')) install.packages('googlesheets4')
if (!require('rlist')) install.packages('rlist')
library(rjson)
library(googlesheets4)
library(rlist)


# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Raw_Data/Precision_Task')
saveDir <- file.path(parentDir,'Data_Compiled/')


setwd(dataDir)

# --- Load Sj Info From Sheets ---
gs4_auth(email = 'jordangarrett@ucsb.edu')
exp_info <- read_sheet('https://docs.google.com/spreadsheets/d/1CtW0BKpcAn0M8aK9PbRPH2Tz7zLYz0nSRQ7KenTdroU/edit#gid=0')


subjects <- exp_info$`Sj ID`
precision_fileNumbs <- exp_info$`File Suffix #...15`

subjects <- subjects[!is.na(subjects)]

exclude_subs <- c(1:4,9:10) #this is key for deciding which subjects are processed

subjects <- setdiff(subjects,exclude_subs)

loc_bins <- c(0,60,120,180,240,300)


# ---------- Custom Functions -------------

"
Adam, Vogel, & Awh, 2017 used model free circular statistics when quatifiying error.
Response error distributions are centered around 0 degrees of error in a circular 
normal space (e.g., -180 == 180). They used functions from 'CircStat' toolbox in Matlab
for stats.
"
# based on the toolbox 'CircStat' from Matlab
compute_circularStats <- function(degrees_error){
  
  stats <- list()
  
  # convert angles to radians
  error.rad <- abs(degrees_error) * pi/180
  
  # circ stats sums exponential of response error in radians times 1i
  error.complex_sum <- sum(exp(1i*error.rad))
  
  # mean is obtained by computing the phase angle of the error vector
  # Arg(z) is the equivalent of angle(z) in Matlab
  # circ_mean.m function in matlab
  error.mean <- Arg(error.complex_sum)
  
  
  # compute mean resultant vector length (MRVL)
  # MRVL indicates variability of the data, 
  # ranging from 0 (complete absence of info about target) to 1 (perfect info about the target)
  # circ_r.m function in matlab
  error.MRVL <- abs(error.complex_sum)/length(error.rad)
  
  stats$mean <- error.mean * 180/pi #convert to degrees
  stats$MRVL <- error.MRVL
  
  return(stats)
}



# ----------- Loop through Sj Files --------

# making this static to ensure the index doesn't change across subjects
all_combos <- list(c(0),c(0,60),c(0,120),c(0,180),c(0,240),c(0,300),
                   c(60),c(60,120),c(60,180),c(60,240),c(60,300),
                   c(120),c(120,180),c(120,240),c(120,300),
                   c(180),c(180,240),c(180,300),
                   c(240),c(240,300),
                   c(300),
                   c(0,60,120,180,240,300))

# check if we have already ran previous subjects. If so, then load that file
if(file.exists(file.path(saveDir,'All_Sj_Precision_Data.RDS'))){
  all_sj.precision_data <- readRDS(file.path(saveDir,'All_Sj_Precision_Data.RDS'))
} else {
  all_sj.precision_data <- list()
}


jatos_file.Prefix <- 'jatos_results_'

for(iSub in 1:length(subjects)){
  
  sj_numb <- subjects[iSub]
  
  # if subject was previously process, skip them
  if(sj_numb %in% all_sj.precision_data$sj_numbs){
    next
  }
  
  # in case we haven't ran the second session for this sj
  if(is.na(precision_fileNumbs[sj_numb])){
    next
  }
  
  file_suffix <- as.character(precision_fileNumbs[sj_numb])
  
  precision_filename <- paste(jatos_file.Prefix,
                              file_suffix,
                              '.txt',sep='')
  
  # load JSON, output is a list of lists
  trial_data <- fromJSON(file = precision_filename)
  
  sample_loc.bins <- trial_data$SampleLocBin
 
  
  #---------------- Overall Precision ----------------------
  sample_orientations <- unlist(trial_data$SampleProbeAngles)
  resp_orientations <- unlist(trial_data$RespProbeAngles)
  
  overall.resp_error <- sample_orientations - resp_orientations
  
  overall.stats <- compute_circularStats(overall.resp_error)
  
  overall_offset_error.mean <- overall.stats$mean
  overall_offset_error.MRVL <- overall.stats$MRVL
  
  # ------------- Precision per location Combo -------------
  # Note: each combo length will tell us the set size, so can use that to parse out set size x combo precision
  
  # there should be an equal number of trials per location combination for comparison (currently =20)
  all_combo.trial_idx <- matrix(nrow = length(all_combos), 
                            ncol = length(trial_data$Trial)/length(all_combos))
  combo.overallError_mean <- c()
  combo.overallError_MRVL <- c()
  for (iLoc_comb in 1:length(all_combos)){
    current_combo <- as.numeric(unlist(all_combos[iLoc_comb]))
    n_combo <- length(list.search(sample_loc.bins, 
                                  identical(.,current_combo)))
    
    all_combo.trial_idx[iLoc_comb,] <- list.findi(sample_loc.bins, 
                                              identical(., current_combo), n=n_combo)
    
    combo.samples <- sample_orientations[all_combo.trial_idx[iLoc_comb,]]
    combo.responses <- resp_orientations[all_combo.trial_idx[iLoc_comb,]]
    
    combo.overall_error <- combo.samples - combo.responses
    
    combo.overall_stats <- compute_circularStats(combo.overall_error)
    
    combo.overallError_mean[iLoc_comb] <- combo.overall_stats$mean
    combo.overallError_MRVL[iLoc_comb] <- combo.overall_stats$MRVL
  }
  
  
  # --------------- Precision per set size ----------------- 
  # Extract trials for each set size
  set_sizes <- sort(unique(trial_data$SetSize))
  
  set_size.error_mean <- c()
  set_size.error_MRVL <- c()
  for(iSet_size in 1:length(set_sizes)){
    
    current.set_size <- set_sizes[iSet_size]
    set_sizeTrials <- which(trial_data$SetSize == current.set_size)
    
    set_size.samples <- unlist(trial_data$SampleProbeAngles[set_sizeTrials])
    set_size.responses <- unlist(trial_data$RespProbeAngles[set_sizeTrials])
    
    set_size.error <- set_size.samples - set_size.responses
    
    set_size.stats <- compute_circularStats(set_size.error)
    
    set_size.error_mean[iSet_size] <- set_size.stats$mean
    set_size.error_MRVL[iSet_size] <- set_size.stats$MRVL
  }
  
  
  # ----------------- Location precision --------------------------
  # Overall precision for each location
  overall_locError.mean <- c()
  overall_locError.MRVL <- c()
  for (iLoc in 1:length(loc_bins)){
    current_loc <- loc_bins[iLoc]
    
    loc.trial_idx <- which(lapply(sample_loc.bins,function(x){is.element(current_loc,x)}) == TRUE)
    
    # grab subset of trials that contain location
    loc_subset_trials <- sample_loc.bins[loc.trial_idx]
    
    # get index of location on subset trials
    loc.trial_pos <- unlist(lapply(sample_loc.bins[loc.trial_idx],function(x){match(current_loc,x)}))
    
    loc.samp_orientations <- trial_data$SampleProbeAngles[loc.trial_idx]
    loc.resp_orientations <- trial_data$RespProbeAngles[loc.trial_idx]
    
    loc.error <- c(rep(0,length(loc_subset_trials)))
    for (iTrial in 1:length(loc_subset_trials)){
      loc.samp_probe <- unlist(loc.samp_orientations[iTrial])[loc.trial_pos[iTrial]]
      loc.resp_probe <- unlist(loc.resp_orientations[iTrial])[loc.trial_pos[iTrial]]
      loc.error[iTrial] <- loc.samp_probe - loc.resp_probe
    }
    
    overall_locStats <- compute_circularStats(loc.error)
    overall_locError.mean[iLoc] <- overall_locStats$mean
    overall_locError.MRVL[iLoc] <- overall_locStats$MRVL
    
  }
  
  # ------------------ Error at location x Set Size -------------------------
  
  # Loop through locations and then each location combo
  # Get error for location of interest when location j present
  ss2.combo_idx <- which(lapply(all_combos,
                                function(x){length(x) == 2}) == TRUE)
  ss2.combos <- all_combos[ss2.combo_idx]
  
  ss2.combo_trials <- all_combo.trial_idx[ss2.combo_idx,]
  
  ss2.raw_error <- data.frame()
  for(iTrial in 1:length(ss2.combo_trials)){
    current_trial <- ss2.combo_trials[iTrial]
    current_combo <- trial_data$SampleLocBin[current_trial]
    sample_orient <- unlist(trial_data$SampleProbeAngles[current_trial])
    resp_orient <- unlist(trial_data$RespProbeAngles[current_trial])
    
    ss2.raw_error[iTrial,c(1:2)] <- sample_orient - resp_orient
    ss2.raw_error[iTrial,3][[1]] <- current_combo
  }
  
  colnames(ss2.raw_error) <- c('Loc1_error','Loc2_error',
                               'LocCombo')
  
  confusion_ss2.mean <- confusion_ss2.MRVL <- matrix(nrow = length(loc_bins), 
                                           ncol = length(loc_bins))
  
  for (iLoc_comb in 1:length(ss2.combos)){
      current_comb <- as.numeric(unlist(ss2.combos[iLoc_comb]))
      comb_idx <- which(lapply(ss2.raw_error[,3], 
                               function(x){identical(x,current_comb)}) == TRUE)
      comb_df <- ss2.raw_error[comb_idx,c(1:2)]
      
      comb_stats <- sapply(comb_df,compute_circularStats)
      
      # index of locations in combo in location bin vector
      matrix_pos <- match(current_comb,loc_bins) 
      
      confusion_ss2.mean[matrix_pos[1],matrix_pos[2]] <- comb_stats[[1]]
      confusion_ss2.mean[matrix_pos[2],matrix_pos[1]] <- comb_stats[[3]]
      
      confusion_ss2.MRVL[matrix_pos[1],matrix_pos[2]] <- comb_stats[[2]]
      confusion_ss2.MRVL[matrix_pos[2],matrix_pos[1]] <- comb_stats[[4]]
      
  }
  
  colnames(confusion_ss2.mean) <- colnames(confusion_ss2.MRVL) <- loc_bins
  rownames(confusion_ss2.mean) <- rownames(confusion_ss2.MRVL) <- loc_bins
  
  
  # Set Size 1 and 6 will just be vectors
  ss1_loc.mean_error <- ss1_loc.MRVL <- c()
  for (iLoc in 1:length(loc_bins)){
    current_loc <- loc_bins[iLoc]
    ss1_loc.trials_idx <- which(lapply(trial_data$SampleLocBin,
                                       function(x){identical(unlist(x),current_loc)}) == TRUE)
    sample_orient <- unlist(trial_data$SampleProbeAngles[ss1_loc.trials_idx])
    resp_orient <- unlist(trial_data$RespProbeAngles[ss1_loc.trials_idx])
    
    error <- sample_orient - resp_orient
    
    ss1_loc.stats <- compute_circularStats(error)
    
    ss1_loc.mean_error[iLoc] <- ss1_loc.stats$mean
    ss1_loc.MRVL[iLoc] <- ss1_loc.stats$MRVL
    
  }
  
  ss6_loc.trials_idx <- which(lapply(trial_data$SampleLocBin,
                              function(x){length(unlist(x)) == 6}) == TRUE)
  ss6_loc.raw_error <- matrix(nrow=length(ss6_loc.trials_idx), ncol = 6)
  for (iTrial in 1:length(ss6_loc.trials_idx)){
    current_trial <- ss6_loc.trials_idx[iTrial]
    
    sample_orient <- unlist(trial_data$SampleProbeAngles[current_trial])
    resp_orient <- unlist(trial_data$RespProbeAngles[current_trial])
    
    ss6_loc.raw_error [iTrial,] <- sample_orient - resp_orient
    
  }
  ss6_loc.stats <- apply(ss6_loc.raw_error, 2,compute_circularStats)
  ss6_loc.mean_error <- unlist(lapply(ss6_loc.stats, function(x){x$mean}))
  ss6_loc.MRVL <- unlist(lapply(ss6_loc.stats, function(x){x$MRVL}))
  
  # ---------------- Store Error Data ----------------
  
  #overall precision
  all_sj.precision_data$overall$mean_error[iSub] <- overall_offset_error.mean
  all_sj.precision_data$overall$MRVL[iSub] <- overall_offset_error.MRVL
  
  #overall precision per combo
  all_sj.precision_data$combo_overall$mean_error[[iSub]] <- combo.overallError_mean
  all_sj.precision_data$combo_overall$MRVL[[iSub]] <- combo.overallError_MRVL
  
  #overall precision per set size
  all_sj.precision_data$ss_overall$mean_error[[iSub]] <- set_size.error_mean
  all_sj.precision_data$ss_overall$MRVL[[iSub]] <- set_size.error_MRVL 
  
  #overall precision per location
  all_sj.precision_data$location_overall$mean_error[[iSub]] <- overall_locError.mean
  all_sj.precision_data$location_overall$MRVL[[iSub]] <- overall_locError.MRVL
  
  #error per set size
  all_sj.precision_data$location.by.ss$ss1$mean_error[[iSub]] <- ss1_loc.mean_error
  all_sj.precision_data$location.by.ss$ss1$MRVL[[iSub]] <- ss1_loc.MRVL
  
  all_sj.precision_data$location.by.ss$ss2$mean_error[[iSub]] <- confusion_ss2.mean
  all_sj.precision_data$location.by.ss$ss2$MRVL[[iSub]] <- confusion_ss2.MRVL
  
  all_sj.precision_data$location.by.ss$ss6$mean_error[[iSub]] <- ss6_loc.mean_error
  all_sj.precision_data$location.by.ss$ss6$MRVL[[iSub]] <- ss6_loc.MRVL
  
  
}

all_sj.precision_data$sj_numbs <- subjects
all_sj.precision_data$loc_combos <- all_combos


# save data 
saveRDS(all_sj.precision_data, file = file.path(saveDir,'All_Sj_Precision_Data.RDS'))



