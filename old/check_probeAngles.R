if (!require('rjson')) install.packages('package')
library(rjson)

dataDir <- '/Users/owner/Downloads/'
result <- fromJSON(file = "jatos_results_20210108003931")

samp.probe_angles <- unlist(result$SampleProbeAngles)
resp.probe_angles <- unlist(result$RespProbeAngles)

error <- abs(samp.probe_angles - resp.probe_angles)

print(paste0('Avg Error: ',round(mean(error),2),'\u00B0'))
print(paste0('Stdv Error: ',round(sd(error),2),'\u00B0'))