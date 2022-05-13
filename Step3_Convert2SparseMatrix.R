#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
#DataSource <- args[1]
DataSource <- 'Dx'

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


tmp <- read.csv(file='/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/Dx_CoOccurMatrix_1.csv')
windows <- as.numeric(names(table(tmp$window)))



## select the max (total) number of code ##

EhrLongFormat <- read.csv('/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/LongFormat_Num/codeRecord.csv')
maxInd <- max(EhrLongFormat$CId)

## merger all 20 chunks ##

for(j in 1:20){
  print(c('chunk: ', j))
  print(Sys.time())
  coccur <- read.csv(paste0('/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/', DataSource, '_CoOccurMatrix_', j, '.csv'))
  
  #if(colnames(coccur)[1]!='code1') {print(j); error_chunk <- c(error_chunk, j); next}
  
  matrices <- list()
  for(i in 1:length(windows)){
    # print(c('window: ', i))
    coccurtemp <- subset(coccur, window==windows[i])
    matrices <- append(matrices, 
                       sparseMatrix(i=as.integer(coccurtemp$code1),
                                    j=as.integer(coccurtemp$code2),
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

save(matrices, windows, file=paste0('/nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/', DataSource, '_CoOccurMatrix.RData'))

