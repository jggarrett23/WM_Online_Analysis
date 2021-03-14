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
    
    Note: hit rate = number of hits / number of change trials
          false alarm rate = number of false alarms/ number of same trials

Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"

# Load necessary packages
if (!require('rjson')) install.packages('rjson')
if (!require('googlesheets4')) install.packages('googlesheets4')
library(rjson)
library(googlesheets4)


# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online/'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Raw_Data/K_Task')
saveDir <- file.path(parentDir,'Data_Compiled')


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
  hits.missing_setSizes <- as.character(setdiff(c(1:6),
                                   as.integer(unlist(dimnames(set_size.hits)))))
  
  
  ## NEED TO FIX HOW WE ARE INPUTTING ZEROS FOR THESE SET SIZES
  
  fas.missing_setSizes <- as.character(setdiff(c(1:6)
                                  ,as.integer(unlist(dimnames(set_size.fas))))) 
  
  set_size.hits <- as.data.frame(set_size.hits)
  set_size.fas <- as.data.frame(set_size.fas)
  
  set_size.hits$Var1 <- factor(set_size.hits$Var1, levels = c(1:6))
  set_size.fas$Var1 <- factor(set_size.fas$Var1, levels = c(1:6))
  
  if (length(hits.missing_setSizes)){
    hits.missing_setSizes <- strsplit(unlist(hits.missing_setSizes),"")
    for(iSet in hits.missing_setSizes){
      set_size.hits <- rbind(set_size.hits, c(iSet,0))
    }
  }
  
  if (length(fas.missing_setSizes)){
    fas.missing_setSizes <- strsplit(unlist(fas.missing_setSizes),"")
    for(iSet in fas.missing_setSizes){
      set_size.fas <- rbind(set_size.fas, c(iSet,0))
    }
  }
  
  set_size.hits <- set_size.hits[order(set_size.hits$Var1),]
  set_size.fas <-  set_size.fas[order(set_size.fas$Var1),]
  
  set_size.hits$Freq <- as.numeric(set_size.hits$Freq) 
  set_size.fas$Freq <- as.numeric(set_size.fas$Freq)
  
  
  hit_rate <- cbind(set_size.hits$Var1,
                                set_size.hits$Freq/n_changeTrials.perSetSize)
  fa_rate <- cbind(set_size.fas$Var1,
                               set_size.fas$Freq/n_sameTrials.perSetSize)
  
  k_df <- as.data.frame(cbind(hit_rate,fa_rate[,2]))
  colnames(k_df) <- c('Set_Size','Hr','FAs')
  
  # N[(h-f)/(1-f)]
  k_df$k_hat <- k_df$Set_Size*((k_df$Hr-k_df$FAs)/(1-k_df$FAs))
  
  colnames(k_df) <- c('Set_Size','Hr','FA','K')
  
  return(k_df)
}

# --- Extract Sj Numbers from Google Sheets ---
gs4_auth(email = 'jordangarrett@ucsb.edu')
exp_info <- read_sheet('https://docs.google.com/spreadsheets/d/1CtW0BKpcAn0M8aK9PbRPH2Tz7zLYz0nSRQ7KenTdroU/edit#gid=0')
k_fileNumbs <- exp_info$`File Suffix #...9`
subjects <- exp_info$`Sj ID`

subjects <- subjects[!is.na(subjects)]


# for sessions that didn't work out (see df notes)
exclude_subs <- c(2)

subjects <- setdiff(subjects,exclude_subs)
# --- Load file containing already calculated Sj Estimates ---
if (file.exists(file.path(dataDir,'All_K_Estimates.csv'))){
  master_k_df <- read.csv(file.path(dataDir,'All_K_Estimates.csv'))
  old_subjects <- master_k_df$Sj_Numb
  
  # only keep subjects that havent been analyzed
  subjects <- setdiff(subjects,old_subjects)
  
} else {
  master_k_df <- data.frame(Sj_Numb = double(),
                            Set6_K = double(),
                            Kmax = double())
}

if(length(subjects) == 0){
  stop('No new subjects to be analyzed')
}
 


# ------------- Extract K Estimates -------------

jatos_file.Prefix <- 'jatos_results_'

# Loop through and extract data from files
all.K_estimates <- data.frame(Sj_Numb = double(),
                              Set6_K = double(),
                              Kmax = double())
for(iSub in 1:length(subjects)){
  
  sj_numb <- subjects[iSub]
  
  # grab file suffix from google sheet
  file_suffix <- as.character(k_fileNumbs[iSub])
  
  k_filename <- paste(jatos_file.Prefix,file_suffix,
                      '.txt',sep = '')
  
  # load JSON, output is a list of lists
  trial_data <- fromJSON(file = k_filename)
  
  # change/no change trial assignments
  trial_type <- trial_data$TrialType
  
  set_sizes <- trial_data$SetSize
  
  # change/no change responses 
  trial_responses <- trial_data$ChangeResponse
  
  # estimate hr, fa, and k
  sdt.k.estimates <- Estimate_K(trial_type,trial_responses,set_sizes)
  
  Kmax <- round(max(sdt.k.estimates['K']),2)
  
  all.K_estimates[iSub,'Sj_Numb'] <- sj_numb
  all.K_estimates[iSub,'Set6_K'] <- round(sdt.k.estimates[6,'K'],2)
  all.K_estimates[iSub,'Kmax'] <- Kmax
}

master_k_df <- rbind(master_k_df,all.K_estimates)

write.csv(master_k_df,
          file.path(saveDir,'All_K_Estimates.csv'),
          row.names = F)

