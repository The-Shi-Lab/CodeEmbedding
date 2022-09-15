#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)


###This is an example of embedding codes with using co-occurrence matrix
###and evaluating the performance different embedding dimentions using
###AUC and NMI


#########################################################
## loading package and function from Xu
#########################################################

library(SNFtool)
library(gtools)
library(stringr)
library(pROC)
library(fastDummies)


source('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/code/EmbeddingFunctions_fromXu.R')


#########################################################
## setting parameters and loading data
#########################################################

#DataSource <- 'Dx'
#DataSource <- 'Px'
#DataSource <- 'Lab'
DataSource <- args[1]
print(DataSource)
tot_EmbedVec <- 500
MAX_ITERS <- 2000
windows<-c(1,2,7,10,14,20,30,40,50,60)-1
iter <- args[2]
threshold <- args[3]

#load(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CoOccurMatrix/', DataSource,'_CoOccurMatrix.RData'))
#load(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/', DataSource,'_CoOccurMatrix/CoOccurMatrix.RData'))
load(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/New_CoOccurMatrix/',threshold,'/CoOccurMatrix.RData'))


#########################################################
## generate embedding
#########################################################

for(i in as.numeric(iter)){
  cocOverall <- as.matrix(matrices[[i]])
  #cocOverall[lower.tri(cocOverall, diag=F)] <- 0
  b <- as(cocOverall, "dgTMatrix")
  coccur <- cbind.data.frame(code1=b@i+1, code2=b@j+1, count=b@x)
  
  # coccur$count <- (coccur$count)
  coccur$code1 <- as.character(coccur$code1)
  coccur$code2 <- as.character(coccur$code2)
  
  singletons <- getSingletonTb(coccur)
  pmi <- construct_pmi(coccur, singletons, my.smooth=0.75)
  sppmi <- construct_sppmi(pmi, k=10)
  
  # set max iteration to 2000
  # set embbeding number to 500 (50 will be the same as first 50)
  rslt <- factor_sppmi(sppmi, dim_size=tot_EmbedVec, iters=MAX_ITERS, remove_empty=T, use_sum=F)
  CodeEmbedVec <- rslt$vecs
  
  print(i)
  print(windows[i])
  # print(head(CodeEmbedVec[,1:4]))
  
  file_name <- paste0("New", '_', windows[i], 'days_', tot_EmbedVec, 'vec')
  #file_name <- paste0(DataSource, '_', windows[i], 'days_', tot_EmbedVec, 'vec')
  
  save(CodeEmbedVec, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeEmbedding/',threshold,'/', file_name, '.RData'))
  write.csv(CodeEmbedVec, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeEmbedding/',threshold,'/', file_name, '.csv'), row.names=F)
}

