#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
path = "/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/"

library(data.table)
library(pROC)
library(dplyr)

DataSource <- args[1]
Window_tmp <- args[2]
Window <- as.numeric(Window_tmp)-1
Version <- args[3]

DataSource<-"New"
# Window = 13
Version<-'ICD9ICD10'
threshold <-5

print(c(DataSource, Window, Version, threshold))

########################
# import crosswalk files
DxCrosswalk <- read.table(paste0(path,'CodeCrosswalk/NewCode_num_source_',threshold,'.txt'), header=T)
DxCrosswalk_freq <- read.table(paste0(path,'UKB_NewCode_freq_20220927_',threshold,'.txt'), header=T)
DxCrosswalk$Freq = DxCrosswalk_freq$Freq[match(DxCrosswalk$NewCode_num, DxCrosswalk_freq$NewCode_num)]

# match icd9 and icd10 codes to phecodes
icd9_phecode <- read.csv('/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/data/phecode_icd9_map_unrolled.csv')
icd10_phecode <- read.csv('/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/data/Phecode_map_v1_2_icd10_beta.csv')

colnames(icd9_phecode) <- c('ICD_dx', 'phecode')
icd9_phecode$phecode_round <- floor(icd9_phecode$phecode)
icd9_phecode <- icd9_phecode[order(icd9_phecode$ICD_dx, icd9_phecode$phecode),]
icd9_phecode <- icd9_phecode[!duplicated(icd9_phecode$ICD_dx),]
icd9_phecode <- icd9_phecode[!is.na(icd9_phecode$phecode),]

icd10_phecode <- icd10_phecode[,c('ICD10','PHECODE')]
colnames(icd10_phecode) <- c('ICD_dx', 'phecode')
icd10_phecode$phecode_round <- floor(icd10_phecode$phecode)
icd10_phecode <- icd10_phecode[order(icd10_phecode$ICD_dx, icd10_phecode$phecode),]
icd10_phecode <- icd10_phecode[!duplicated(icd10_phecode$ICD_dx),]
icd10_phecode <- icd10_phecode[!is.na(icd10_phecode$phecode),] ## there are 105 codes with NA phecode

icd_phecode <- rbind(icd9_phecode, icd10_phecode)

############################
# import the embedding matrix
file_name <- paste0(DataSource, '_', Window, 'days_500vec.RData')
load(paste0(path,'CodeEmbedding/', file_name)) # load embedding matrix "CodeEmbedVec"

# determine codes for evaluation
all_codes <- data.frame("Newcode_num"=rownames(CodeEmbedVec))
all_codes1 <- all_codes %>%
  mutate("Source" = DxCrosswalk$Source[match(Newcode_num,as.numeric(DxCrosswalk$NewCode_num)-100000000)]) %>%
  filter(Source=="Original") %>%
  mutate("Newcode" = DxCrosswalk$NewCode[match(Newcode_num,as.numeric(DxCrosswalk$NewCode_num)-100000000)]) %>%
  mutate("Freq" = DxCrosswalk$Freq[match(Newcode_num,as.numeric(DxCrosswalk$NewCode_num)-100000000)]) %>%
  mutate("Newcode1" = sub('Org_', '', Newcode)) %>%
  mutate("Phecode_round" = icd_phecode$phecode_round[match(Newcode1, icd_phecode$ICD_dx)]) %>%
  filter(Freq>=5)
codes_for_evaluation = na.omit(all_codes1)
dim(codes_for_evaluation) 

# compute the reference vector
print(Sys.time())
ref_vec <- c()
for(i in 2:nrow(codes_for_evaluation)){
  tmp_phecode <- codes_for_evaluation$Phecode_round[i]
  tmp_ref_vec <- (codes_for_evaluation$Phecode_round[1:(i-1)]==tmp_phecode)*1
  
  ref_vec <- c(ref_vec, tmp_ref_vec)
  if(i%%1000==0) print(i)
}
print('ref_vec Done')
print(Sys.time())

# evaluate the prediction performance of the embedding matrix
CodeEmbedVec_evaluate <- CodeEmbedVec[codes_for_evaluation$Newcode_num,] # select rows in the Embedding Vector
print(start)

for(j in c(10,30,50,100,150,200,250,300,350,400,450,500)){
  tmp <- CodeEmbedVec_evaluate[,1:j]
  tmp2 <- sqrt(apply(tmp^2,1,sum))
  tmp_norm <- tmp/tmp2
  cosine_matrix <- tmp_norm%*%t(tmp_norm)
  cosine_vec <- cosine_matrix[upper.tri(cosine_matrix)]
  ROC_result <- roc(ref_vec, cosine_vec)
  
  out_file_name <- paste0('AUC_', DataSource, '_', Window, 'days_', j, 'vec_', Version,'.RData')
  save(ROC_result, file=paste0(path,'Evaluation/',out_file_name))
  print(j)
  print(Sys.time())
}


