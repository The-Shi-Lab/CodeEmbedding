#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
DataSource <- args[1]
threshold <- args[2]
################################################
#                   Load data                  #
################################################

library(dplyr)
library(SNFtool)
library(gtools)
library(stringr)
library(pROC)
library(fastDummies)
library(Matrix)

DxCrosswalk <- read.table(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeCrosswalk/NewCode_num_',threshold,'.txt'), header=T)
##DxCrosswalk <- read.table('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeCrosswalk/DxCode_num.txt', header=T)
#PxCrosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/PxCode_num.txt', header=T)
#LabCrosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/LabCode_num.txt', header=T)

tmp <- read.csv(file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/New_CoOccurMatrix/',threshold,'/CoOccurMatrix_1000.csv'))
#tmp <- read.csv(file='/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Dx_CoOccurMatrix/CoOccurMatrix_1000.csv')
windows <- as.numeric(names(table(tmp$window)))


## select the max (total) number of code ##

if(DataSource=='Dx'){
  maxInd <- nrow(DxCrosswalk)
  RemovePrefix <- 100000000
}
if(DataSource=='Px'){
  maxInd <- nrow(PxCrosswalk)
  RemovePrefix <- 200000000
}
if(DataSource=='Lab'){
  maxInd <- nrow(LabCrosswalk)
  RemovePrefix <- 300000000
}



## merge all 3000 chunks ##

for(j in 1:3000){
  print(c('chunk: ', j))
  print(Sys.time())
  coccur <- read.csv(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/New_CoOccurMatrix/',threshold,'/CoOccurMatrix_', j, '.csv'))
  
  #if(colnames(coccur)[1]!='code1') {print(j); error_chunk <- c(error_chunk, j); next}
  
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

save(matrices, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/New_CoOccurMatrix/',threshold,'/CoOccurMatrix.RData'))
# save(matrices, windows, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CoOccurMatrix/Px_CoOccurMatrix.RData')
# save(matrices, windows, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CoOccurMatrix/Lab_CoOccurMatrix.RData')
#load("/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/New_CoOccurMatrix/CoOccurMatrix.RData")
