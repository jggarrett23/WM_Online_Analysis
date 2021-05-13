"
-----------------------------------------------------------
Compile all Sj error estimates


Author: Jordan Garrett
UCSB Attention Lab
jordangarrett@ucsb.edu
-----------------------------------------------------------
"
# Load necessary packages
library(readr)

# ------------- Set up Directories -------------
if (Sys.info()['sysname'] == 'Windows'){
  parentDir <- 'D:/WM_Online'
} else{
  parentDir <- '/Users/owner/Downloads'
}

dataDir <- file.path(parentDir,'Data/Precision/')
saveDir <- file.path(parentDir,'Data_Compiled/')


setwd(dataDir)


# ---------- Grab Sj Files ---------------

sj_files <- list.files(dataDir)
subjects <- unlist(lapply(sj_files,parse_number))


# --------- Loop through and Compile --------
all_sj.precision_data <- c()

for (iSub in 1:length(subjects)){
  sj_numb <- subjects[iSub]
  
  sj.file_name <- sj_files[iSub]
  
  sj.data <- readRDS(file.path(dataDir,
                               sprintf('Sj%02d_Precision_Data.RDS',sj_numb)))
  
  #overall precision
  all_sj.precision_data$overall$mean_error[iSub] <- sj.data$overall$mean_error
  all_sj.precision_data$overall$MRVL[iSub] <- sj.data$overall$MRVL
  
  #overall precision per combo
  all_sj.precision_data$combo_overall$mean_error[[iSub]] <- sj.data$combo_overall$mean_error
  all_sj.precision_data$combo_overall$MRVL[[iSub]] <- sj.data$combo_overall$MRVL
  
  #overall precision per set size
  all_sj.precision_data$ss_overall$mean_error[[iSub]] <- sj.data$ss_overall$mean_error
  all_sj.precision_data$ss_overall$MRVL[[iSub]] <- sj.data$ss_overall$MRVL
  
  #overall precision per location
  all_sj.precision_data$location_overall$mean_error[[iSub]] <- sj.data$location_overall$mean_error
  all_sj.precision_data$location_overall$MRVL[[iSub]] <- sj.data$location_overall$MRVL
  
  #error per location x set size
  all_sj.precision_data$location.by.ss$ss1$mean_error[[iSub]] <- sj.data$location.by.ss$ss1$mean_error
  all_sj.precision_data$location.by.ss$ss1$MRVL[[iSub]] <- sj.data$location.by.ss$ss1$MRVL
  
  all_sj.precision_data$location.by.ss$ss2$mean_error[[iSub]] <- sj.data$location.by.ss$ss2$mean_error
  all_sj.precision_data$location.by.ss$ss2$MRVL[[iSub]] <- sj.data$location.by.ss$ss2$MRVL
  
  all_sj.precision_data$location.by.ss$ss6$mean_error[[iSub]] <- sj.data$location.by.ss$ss6$mean_error
  all_sj.precision_data$location.by.ss$ss6$MRVL[[iSub]] <- sj.data$location.by.ss$ss6$MRVL
  
  all_sj.precision_data$location.by.ss$overall$mean_error[[iSub]] <- sj.data$location.by.ss$overall$mean_error
  all_sj.precision_data$location.by.ss$overall$MRVL[[iSub]] <- sj.data$location.by.ss$overall$MRVL
  
  #categorical confusion matrices
  all_sj.precision_data$cc_matrix[[iSub]] <- sj.data$cc_matrix
}

all_sj.precision_data$sj_numbs <- subjects
all_sj.precision_data$loc_combos <- sj.data$loc_combos

# save data 
saveRDS(all_sj.precision_data, 
        file = file.path(saveDir,'All_Sj_Precision_Data.RDS'))


