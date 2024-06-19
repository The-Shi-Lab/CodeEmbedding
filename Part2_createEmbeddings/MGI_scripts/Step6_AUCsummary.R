#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(pROC)
library(reshape)
library(ggplot2)
library(plotly) 
library(hrbrthemes)

DataSource <- args[1]
i <- args[2]
Version <- args[3]
DataSource <- 'New'
Version <-"ICD9ICD10"

path = "/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/Evaluation/"
# Windows <- c(13,19,29,39,49,59,69,79,89,99)
Windows <- c(0, 1, 6, 9, 13, 19, 29, 39, 49, 59)
print(c(DataSource, Version, i))

out <- c()
for(j in Windows){
  file_name <- paste0('AUC_', DataSource, '_', j, 'days_', i, 'vec_', Version, '.RData')
  file.exists(paste0(path, file_name))
  load(paste0(path, file_name))

  tmp <- ci(ROC_result)
  out <- rbind(out, c(tmp[2], tmp[1], tmp[3]))
  print(c(j, 'daysrbind(out, c(tmp[2], t'))
}
out_name <- paste0('AUC_', DataSource, i, 'vec_', Version, '.csv')
out_name
path
write.csv(out, file=paste0(path,'NewAUCsummary/', out_name),row.names=F) 



