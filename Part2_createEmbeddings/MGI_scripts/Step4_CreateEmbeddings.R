#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)


if("SNFtool" %in% rownames(installed.packages()) == FALSE) {install.packages("SNFtool")}
if("gtools" %in% rownames(installed.packages()) == FALSE) {install.packages("gtools")}
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}
if("pROC" %in% rownames(installed.packages()) == FALSE) {install.packages("pROC")}
if("fastDummies" %in% rownames(installed.packages()) == FALSE) {install.packages("fastDummies")}

library(SNFtool)
library(gtools)
library(stringr)
library(pROC)
library(fastDummies)

# loading package and embedding-creating function from Xu
source('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/scripts/EmbeddingFunctions_fromXu.R')

# set parameters and load the sparse cooccur matrix
DataSource <- 'New'
threshold <- 10

DataSource <- args[1]
iter <- args[2]

print(DataSource)
tot_EmbedVec <- 500
MAX_ITERS <- 2000
# windows<-c(1,2,7,10,14,20,30,40,50,60)-1
windows<-c(0,1,6,9,13,19,29,39,49,59)

path = "/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/"

load(paste0(path,'New_CoOccurMatrix/CoOccurMatrix.RData'))

# generate embedding
set.seed(2022)

for(i in as.numeric(iter)){
  cocOverall <- as.matrix(matrices[[i]])
  b <- as(cocOverall, "dgTMatrix")
  coccur <- cbind.data.frame(code1=b@i+1, code2=b@j+1, count=b@x)
  
  coccur$code1 <- as.character(coccur$code1)
  coccur$code2 <- as.character(coccur$code2)
  
  singletons <- getSingletonTb(coccur)
  pmi <- construct_pmi(coccur, singletons, my.smooth=0.75)
  sppmi <- construct_sppmi_pmi(pmi, k=10)
  
  # set max iteration to 2000
  # set embbeding number to 500 (50 will be the same as first 50)
  rslt <- factor_sppmi(sppmi, dim_size=tot_EmbedVec, iters=MAX_ITERS, remove_empty=T, use_sum=F)
  CodeEmbedVec <- rslt$vecs

  print(i)
  print(windows[i])

  file_name <- paste0("New", '_', windows[i], 'days_', tot_EmbedVec, 'vec')

  save(CodeEmbedVec, file=paste0(path,'CodeEmbedding/', file_name, '.RData'))
  write.csv(CodeEmbedVec, file=paste0(path,'CodeEmbedding/', file_name, '.csv'), row.names=F)
}

