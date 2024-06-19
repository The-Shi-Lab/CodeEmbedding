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


source('/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/scripts/EmbeddingFunctions_fromXu.R')


#########################################################
## setting parameters and loading data
#########################################################
# 
# DataSource <- 'Dx'
#DataSource <- 'Px'
#DataSource <- 'Lab'
# iter = 1:10
threshold <- 5

DataSource <- args[1]
iter <- args[2]
# threshold <- args[3]
print(DataSource)
tot_EmbedVec <- 500
MAX_ITERS <- 2000
# windows<-c(1,2,7,10,14,20,30,40,50,60)-1
# windows<-c(13,19,29,39,49,59,69,79,89,99)
windows<-c(69,79,89,99,109,119,129,139,149,159)

path = "/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/"

load(paste0(path,'New_CoOccurMatrix/CoOccurMatrix.RData'))


#########################################################
## generate embedding
#########################################################

set.seed(2022)

# "iter" represents the window #
for(i in as.numeric(iter)){
  cocOverall <- as.matrix(matrices[[i]])
  b <- as(cocOverall, "dgTMatrix")
  coccur <- cbind.data.frame(code1=b@i+1, code2=b@j+1, count=b@x)
  
  coccur$code1 <- as.character(coccur$code1)
  coccur$code2 <- as.character(coccur$code2)
  
  singletons <- getSingletonTb(coccur)
  pmi <- construct_pmi(coccur, singletons, my.smooth=0.75)
  sppmi <- construct_sppmi(pmi, k=10) # check
  
  # set max iteration to 2000
  # set embbeding number to 500 (50 will be the same as first 50)
  rslt <- factor_sppmi(sppmi, dim_size=tot_EmbedVec, iters=MAX_ITERS, remove_empty=T, use_sum=F)
  CodeEmbedVec <- rslt$vecs # check
  
  print(i)
  print(windows[i])

  file_name <- paste0("New", '_', windows[i], 'days_', tot_EmbedVec, 'vec')

  save(CodeEmbedVec, file=paste0(path,'CodeEmbedding/', file_name, '.RData'))
  write.csv(CodeEmbedVec, file=paste0(path,'CodeEmbedding/', file_name, '.csv'), row.names=F)
}

