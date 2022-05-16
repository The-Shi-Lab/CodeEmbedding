#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(pROC)

DataSource <- args[1]
Window_tmp <- args[2]
Window <- as.numeric(Window_tmp)-1
Version <- args[3]

print(c(DataSource, Window, Version))

##############
### Part 1 ###
##############

DxCrosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding_SurgPcpCohort/data/CodeCrosswalk/DxCode_num.txt', header=T)
dim(DxCrosswalk)
# [1] 27014     2

icd9_phecode <- read.csv('/nfs/turbo/mgi-shixu/data/MGI/dx/phecode_icd9_map_unrolled.csv')
icd10_phecode <- read.csv('/nfs/turbo/mgi-shixu/data/MGI/dx/Phecode_map_v1_2_icd10_beta.csv')

colnames(icd9_phecode) <- c('ICD_dx', 'phecode')
icd9_phecode$phecode_round <- round(icd9_phecode$phecode)
icd9_phecode <- icd9_phecode[order(icd9_phecode$ICD_dx, icd9_phecode$phecode),]
icd9_phecode <- icd9_phecode[!duplicated(icd9_phecode$ICD_dx),]

icd10_phecode <- icd10_phecode[,c('ICD10','PHECODE')]
colnames(icd10_phecode) <- c('ICD_dx', 'phecode')
icd10_phecode$phecode_round <- round(icd10_phecode$phecode)
icd10_phecode <- icd10_phecode[order(icd10_phecode$ICD_dx, icd10_phecode$phecode),]
icd10_phecode <- icd10_phecode[!duplicated(icd10_phecode$ICD_dx),]

if(Version=='ICD9'){icd_phecode <- icd9_phecode}
if(Version=='ICD10'){icd_phecode <- icd10_phecode}
if(Version=='ICD9ICD10'){icd_phecode <- rbind(icd9_phecode, icd10_phecode)}


ICDinDx <- icd_phecode[icd_phecode$ICD_dx%in%DxCrosswalk$DxCode,]
idx <- match(ICDinDx$ICD_dx, DxCrosswalk$DxCode)
ICDinDx$DxCode_num <- DxCrosswalk$DxCode_num[idx] - 100000000
# 7412 ICD-9
dim(ICDinDx)
# [1] 10656     4

# in SurgPcpCohort #
#ICD9: 7708
#ICD10: 3424
#ICD9ICD10: 11132



##############
### Part 2 ###
##############

# DataSource <- 'Dx'
# Window <- 7
file_name <- paste0(DataSource, '_', Window, 'days_500vec.RData')
load(paste0('/nfs/turbo/mgi-shixu/project/CodeEmbedding_SurgPcpCohort/data/CodeEmbedding/', file_name))
dim(CodeEmbedVec)
# [1] 26964   500

CodeInCE <- rownames(CodeEmbedVec)

VecInPhe <- as.numeric(CodeInCE)%in%ICDinDx$DxCode_num
table(VecInPhe)
# VecInPhe
# FALSE  TRUE 
# 16346 10618 
CodeEmbedVec_tmp <- CodeEmbedVec[VecInPhe,]

ICDinDx_formatted <- ICDinDx[match(rownames(CodeEmbedVec_tmp), ICDinDx$DxCode_num),]
dim(ICDinDx_formatted)
# [1] 10618     4
table(rownames(CodeEmbedVec_tmp)==ICDinDx_formatted$DxCode_num)
# TRUE 
# 10618 


print(Sys.time())
SamePhecode <- c()
for(i in 2:nrow(CodeEmbedVec_tmp)){
  tmp_phecode <- ICDinDx_formatted$phecode_round[i]
  tmp_SamePhecode <- (ICDinDx_formatted$phecode_round[1:(i-1)]==tmp_phecode)*1
  
  SamePhecode <- c(SamePhecode, tmp_SamePhecode)
  if(i%%1000==0) print(i)
}
print('SamePhecode Done')
print(Sys.time())

for(j in c(10,30,50,100,150,200,250,300,350,400,450,500)){
  tmp <- CodeEmbedVec_tmp[,1:j]
  tmp2 <- sqrt(apply(tmp^2,1,sum))
  tmp_norm <- tmp/tmp2
  cosine_matrix <- tmp_norm%*%t(tmp_norm)
  cosine_vec <- cosine_matrix[upper.tri(cosine_matrix)]
  
  # COSINE_MATRIX <- cbind(COSINE_MATRIX, cosine_vec)
  ROC_result <- roc(SamePhecode, cosine_vec)
  # AUC_result <- auc(ROC_result)
  # assign(paste0('AUC_', j, 'vec'), AUC_result)
  
  out_file_name <- paste0('AUC_', DataSource, '_', Window, 'days_', j, 'vec_', Version,'.RData')
  # save(ROC_result, AUC_result, file=paste0('/nfs/turbo/mgi-shixu/project/CodeEmbedding/result/CosinePhecodeAUC/', out_file_name))
  save(ROC_result, file=paste0('/nfs/turbo/mgi-shixu/project/CodeEmbedding_SurgPcpCohort/result/CosinePhecodeAUC/', out_file_name))

  print(j)
  print(Sys.time())
}


