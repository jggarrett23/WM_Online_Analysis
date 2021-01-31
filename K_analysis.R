"
-----------------------------------------------------------
Extract K estimates from whole report change detection task.
K algorithm used from Pashler (1988):
  N[(h-f)/(1-f)]
  
  where 
    N = set size
    h = hit rate (Hit = correct change detection; 
                  Correct Rejection = correct no change detection)
    f = false alarm rate (Responding change when no change)

Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"

# Load necessary packages
if (!require('rjson')) install.packages('rjson')
library(rjson)


# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online/'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'K_Task_Results')


setwd(dataDir)

# ------------- Create function to estimate K -------------
Estimate_K <- function(trial_type,responses,set_sizes){
  change_trials <- which(trial_type == 1)
  same_trials <- which(trial_type == 0)
  
  change_responseTrials <- which(responses == 1)
  
  # intersect will give us overlapping elements between both vectors,
  # even if they are of different length
  hitTrials_idx <- intersect(change_trials,change_responseTrials)
  faTrials_idx <- intersect(same_trials,change_responseTrials)
  
  # use the index of hit/fa trials create above and determine the set 
  # size on those trials
  
  # hit rate = # hits / # change trials
  # fa rate = # fa / # same trials
  n_changeTrials.perSetSize <- table(set_sizes[change_trials]) # should be 30 for each
  n_sameTrials.perSetSize <- table(set_sizes[same_trials])
  
  set_size.hits <- table(set_sizes[hitTrials_idx])
  set_size.fas <- table(set_sizes[faTrials_idx])
  
  # check if there are set sizes that need to be logged as zero
  hits.missing_setSizes <- setdiff(dimnames(set_size.hits), 
                                   list(as.character(c(1:6))))
  
  ## NEED TO FIX HOW WE ARE INPUTTING ZEROS FOR THESE SET SIZES
  set_size.hits[array(as.character(hits.missing_setSizes))] <- 0
  
  fas.missing_setSizes <- setdiff(dimnames(set_size.fas), 
                                   list(as.character(c(1:6))))
  
  set_size.fas[as.character(fas.missing_setSizes)] <- 0
  
  set_size.hits <- set_size.hits[order(names(set_size.hits))]
  set_size.fas <- set_size.fas[order(names(set_size.fas))]
  
  hit_df <- as.data.frame(set_size.hits/n_changeTrials.perSetSize)
  fa_df <- as.data.frame(set_size.fas/n_sameTrials.perSetSize)
  
  k_df <- cbind(hit_df,fa_df[2])
  colnames(k_df) <- c('Set_Size','Hr','FAs')
  
  k_df$Set_Size <- as.numeric(k_df$Set_Size)
  
  # N[(h-f)/(1-f)]
  k_df$k_hat <- k_df$Set_Size*((k_df$Hr-k_df$FAs)/(1-k_df$FAs))
  
  colnames(k_df) <- c('Set_Size','Hr','FA','K')
  
  return(k_df)
}

# ------------- Extract K Estimates -------------

# Dont have Sj numbs, so need to extract all the files
k_files <- list.files(dataDir, pattern = "*jordan.txt")

# Loop through and extract data from files
all.K_estimates <- array(dim = length(k_files))
for(iFile in 1:length(k_files)){
  
  # load JSON, output is a list of lists
  trial_data <- fromJSON(file = k_files[iFile])
  
  # change/no change trial assignments
  trial_type <- trial_data$TrialType
  
  set_sizes <- trial_data$SetSize
  
  # change/no change responses 
  trial_responses <- trial_data$ChangeResponse
  
  # estimate hr, fa, and k
  sdt.k.estimates <- Estimate_K(trial_type,trial_responses,set_sizes)
  
  average.k <- round(mean(sdt.k.estimates$K),2)
  
  all.K_estimates[iFile] <- average.k
}
