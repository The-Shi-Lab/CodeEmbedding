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
i <- args[2]
#vecleng= c(10,30,50,100,150,200,250,300,350,400,450,500)
Version <- args[3]
#Version <-"ICD9ICD10"


path = "/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/Evaluation/"
# Windows <- c(13,19,29,39,49,59,69,79,89,99)
Windows<-c(69,79,89,99,109,119,129,139,149,159)

print(c(DataSource, Version, i))


out <- c()
for(j in Windows){
  file_name <- paste0('AUC_', DataSource, '_', j, 'days_', i, 'vec_', Version, '.RData')
  load(paste0(path, file_name))
  tmp <- ci(ROC_result)
  out <- rbind(out, c(tmp[2], tmp[1], tmp[3]))
  print(c(j, 'daysrbind(out, c(tmp[2], t'))
}
out_name <- paste0('AUC_', DataSource, i, 'vec_', Version, '.csv')
write.csv(out, file=paste0(path,'NewAUCsummary/', out_name),row.names=F) 



