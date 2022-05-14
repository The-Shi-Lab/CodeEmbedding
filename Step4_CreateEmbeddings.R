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


source('/nfs/turbo/mgi-shixu/project/CodeEmbedding/code/EmbeddingFunctions_fromXu.R')


#########################################################
## setting parameters and loading data
#########################################################

DataSource <- 'Dx'
#DataSource <- args[1]
print(DataSource)
tot_EmbedVec <- 500
MAX_ITERS <- 2000

load(paste0('/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/', DataSource,'_CoOccurMatrix.RData'))

tmp <- read.csv(file='/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/Dx_CoOccurMatrix_1.csv')
windows <- as.numeric(names(table(tmp$window)))

#iter <- args[2]


#########################################################
## generate embedding
#########################################################

#for(i in as.numeric(iter)){
for(i in 1:length(windows)){
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
  
  file_name <- paste0(DataSource, '_', windows[i], 'days_', tot_EmbedVec, 'vec')
  save(CodeEmbedVec, file=paste0('/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CodeEmbedding/', file_name, '.RData'))
  write.csv(CodeEmbedVec, file=paste0('/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CodeEmbedding/', file_name, '.csv'), row.names=F)
}

