#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(pROC)
library("reshape")
library("ggplot2")
library("plotly") 
library(hrbrthemes)

DataSource <- args[1]
#DataSource <- 'Dx'
#DataSource <- 'New'
i <- args[2]
#vecleng= c(10,30,50,100,150,200,250,300,350,400,450,500)
Version <- args[3]
#Version <-"ICD9ICD10"
threshold <-args[4]

file_path <- paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/NewCosinePhecodeAUC/',threshold, '/')
#file_path <- '/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/DxCosinePhecodeAUC/'
# Windows <- c(1, 7, 10, 14, 20, 30, 40, 50, 60)
Windows <- c(0, 1, 6, 9, 13, 19, 29, 39, 49, 59)

print(c(DataSource, Version, i))


out <- c()
for(j in Windows){
  file_name <- paste0('AUC_', DataSource, '_', j, 'days_', i, 'vec_', Version, '.RData')
  load(paste0(file_path, file_name))
  tmp <- ci(ROC_result)
  out <- rbind(out, c(tmp[2], tmp[1], tmp[3]))
  print(c(j, 'daysrbind(out, c(tmp[2], t'))
}
out_name <- paste0('AUC_', DataSource, i, 'vec_', Version, '.csv')
write.csv(out, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/NewAUCsummary/', threshold,'/', out_name),
          row.names=F) 
#write.csv(out, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/DxAUCsummary/', out_name),row.names=F) 


