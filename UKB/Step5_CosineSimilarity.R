#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(pROC)

DataSource <- args[1]
Window_tmp <- args[2]
Window <- as.numeric(Window_tmp)-1
Version <- args[3]
threshold <- args[4]

#DataSource<-"New"
#Window<-c(1,2,7,10,14,20,30,40,50,60)-1
#Window=0
#Version<-'ICD9ICD10'

print(c(DataSource, Window, Version, threshold))

### Part 1 ###

DxCrosswalk <- read.table(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeCrosswalk/NewCode_num_source_',threshold,'.txt'), header=T)
#DxCrosswalk <- read.table('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeCrosswalk/DxCode_num.txt', header=T)

dim(DxCrosswalk)
# [1] 12624     3 --- method 2
# [1] 15204     2

#####
#dx_freq<-read.table('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/UKB_NewCode_freq_20220708.txt', header=T)
#dx_freq<-read.table('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/UKB_DxCode_freq_20220708.txt', header=T)
#####

icd9_phecode <- read.csv('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/phecode_icd9_map_unrolled.csv')
icd10_phecode <- read.csv('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/Phecode_map_v1_2_icd10_beta.csv')

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

if(Version=='ICD9'){icd_phecode <- icd9_phecode}
if(Version=='ICD10'){icd_phecode <- icd10_phecode}
if(Version=='ICD9ICD10'){icd_phecode <- rbind(icd9_phecode, icd10_phecode)}


ICDinDx <- icd_phecode[icd_phecode$ICD_dx%in%DxCrosswalk$NewCode[DxCrosswalk$Source!="Phecode"],]
dim(ICDinDx)
## [1] 7212    3

###########
#idx2 <- match(ICDinDx$DxCode_num, dx_freq$Dxcode)
#ICDinDx$DxCode_freq <- dx_freq$Freq[idx2]
#ICDinDx$DxCode_num <- ICDinDx$DxCode_num - 10000000
#Phecode_round_Freq=as.data.frame(table(ICDinDx$phecode_round))
#colnames(Phecode_round_Freq)=c("phecode_round","phecode_round_freq")
#idx3 <- match(ICDinDx$phecode_round, Phecode_round_Freq$phecode_round)
#ICDinDx$phecode_round_freq <- Phecode_round_Freq$phecode_round_freq[idx3]
###########


idx <- match(ICDinDx$ICD_dx, DxCrosswalk$NewCode)
ICDinDx$NewCode_num <- DxCrosswalk$NewCode_num[idx]-100000000
ICDinDx$Source<- DxCrosswalk$Source[idx]

dim(ICDinDx)
# [1] 7212    5 --- method 2 grouped up
# [1] 9340    6
### Part 2 ###

# Window <- 0
file_name <- paste0(DataSource, '_', Window, 'days_500vec.RData')
#file_name <- paste0(DataSource, '_', Window, 'days_500vec.RData')
load(paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/CodeEmbedding/',threshold,'/', file_name))
dim(CodeEmbedVec)
# [1] 12129   500 --- method 2
# [1] 14075   500


CodeInCE <- rownames(CodeEmbedVec)

VecInPhe <- (as.numeric(CodeInCE)%in%ICDinDx$NewCode_num[ICDinDx$Source=="Original"])
table(VecInPhe)
# VecInPhe
# FALSE  TRUE 
#  5160  6969 --- method 2
#  5681  8488 --- window6

CodeEmbedVec_tmp <- CodeEmbedVec[VecInPhe,]
########################
# sort our target codes with frequency>=3 similar to what we did towards training codes
#NewCode_freq<-fread("./UKB_icdcode/UKB_NewCode_freq_20220708.txt",sep="\t")
#NewCode_Frequent=NewCode_freq[NewCode_freq$Freq>=3,]
#ICDinDx_Frequent=ICDinDx[ICDinDx$DxCode_freq>=3|ICDinDx$phecode_round_freq>=10,]
#ICDinDx_rare=ICDinDx[!(ICDinDx$DxCode_freq>=3|ICDinDx$phecode_round_freq>=10),]

#VecLargeFreq <- as.numeric(rownames(CodeEmbedVec_tmp)) %in% ICDinDx_Frequent$DxCode_num
#table(VecLargeFreq)
#CodeEmbedVec_tmp <-CodeEmbedVec_tmp[VecLargeFreq,]
########################

ICDinDx_formatted <- ICDinDx[match(rownames(CodeEmbedVec_tmp), ICDinDx$NewCode_num),]
dim(ICDinDx_formatted)
# [1] 6969    5 --- method 2
# [1] 8553    4 --- window0

table(rownames(CodeEmbedVec_tmp)==ICDinDx_formatted$NewCode_num)
# TRUE 
# 6969
# 8553

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
  
  save(ROC_result, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/NewCosinePhecodeAUC/',threshold, '/',
                               out_file_name))
  #save(ROC_result, file=paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/DxCosinePhecodeAUC/', out_file_name))
  
  print(j)
  print(Sys.time())
}


