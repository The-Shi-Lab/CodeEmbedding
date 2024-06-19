#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
DataSource <- args[1]
DataSource = 'New'
threshold=10
# path = "/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10/"
path = "/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/"


if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("SNFtool" %in% rownames(installed.packages()) == FALSE) {install.packages("SNFtool")}
if("gtools" %in% rownames(installed.packages()) == FALSE) {install.packages("gtools")}
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}
if("pROC" %in% rownames(installed.packages()) == FALSE) {install.packages("pROC")}
if("fastDummies" %in% rownames(installed.packages()) == FALSE) {install.packages("fastDummies")}
if("Matrix" %in% rownames(installed.packages()) == FALSE) {install.packages("Matrix")}

library(dplyr)
library(SNFtool)
library(gtools)
library(stringr)
library(pROC)
library(fastDummies)
library(Matrix)

DxCrosswalk <- read.table(paste0(path,'CodeCrosswalk/NewCode_num_',threshold,'.txt'), header=T)
tmp <- read.csv(file=paste0(path,'New_CoOccurMatrix/CoOccurMatrix_1.csv'))

windows <- as.numeric(names(table(tmp$window)))
maxInd <- nrow(DxCrosswalk)
RemovePrefix <- 100000000

## merge all 3000 chunks ##
for(j in 1:3000){
  print(c('chunk: ', j))
  print(Sys.time())
  coccur <- read.csv(paste0(path,'New_CoOccurMatrix/CoOccurMatrix_', j, '.csv'))
  
  if(colnames(coccur)[1]!='code1') {print(j); error_chunk <- c(error_chunk, j); next}
  
  matrices <- list()
  
  for(i in 1:length(windows)){
    # print(c('window: ', i))
    coccurtemp <- subset(coccur, window==windows[i])
    matrices <- append(matrices, 
                       sparseMatrix(i=as.integer(coccurtemp$code1-RemovePrefix),
                                    j=as.integer(coccurtemp$code2-RemovePrefix),
                                    x=coccurtemp$count,
                                    dims=rep(maxInd,2)))
  }
  
  if(j==1){
    matrices_all <- matrices
  }else{
    for(k in 1:length(windows)){
      # print(k)
      matrices_all[[k]] <- matrices_all[[k]]+matrices[[k]]
    }
  }
  
}


## sum up counts for different windows ##
matrices_all_accum <- matrices_all
for(i in 2:length(windows)){
  matrices_all_accum[[i]] <- matrices_all_accum[[i]]+matrices_all_accum[[i-1]]
}

## output result ##
matrices <- matrices_all_accum
save(matrices, file=paste0(path,'New_CoOccurMatrix/CoOccurMatrix.RData'))
