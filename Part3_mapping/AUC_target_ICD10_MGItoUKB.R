library(tidyverse)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(pROC)

## Source datasets
MGI_freq_all_noDuplicates = read.csv("~/EHRmapping/mapping/data/freq_MGI_id.csv")[,-1]
UKB_freq_all_noDuplicates = read.csv("~/EHRmapping/mapping/data/freq_UKB_id.csv")[,-1]
MGI_freq_all_withDuplicates = read.csv("~/EHRmapping/mapping/data/freq_MGI_withDuplicates.csv")[,-1]
UKB_freq_all_withDuplicates = read.csv("~/EHRmapping/mapping/data/freq_UKB_withDuplicates.csv")[,-1]

load("~/EHRmapping/mapping/data/UKB_ICD10_noDuplicates.RData")
load("~/EHRmapping/mapping/data/MGI_ICD10_noDuplicates.RData")
dim(UKB_ICD10)
dim(MGI_ICD10)

source("~/EHRmapping/mapping/scripts/function_mapping.R")
source("~/EHRmapping/mapping/scripts/refine_mapping.R")

#### 
embed_length = 500

######UKB embedding #####
# read New_num for embeddings
load("/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeEmbedding/New_129days_500vec.RData")
CodeEmbedVec_UKB = CodeEmbedVec
row_idx = as.numeric(row.names(CodeEmbedVec_UKB))
# read embeddings
UKB_New_500vec = read_csv("/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeEmbedding/New_129days_500vec.csv")
UKB_New_500vec = cbind(row_idx,UKB_New_500vec)
# read crosswalk file for mapping
crosswalk = read.table('/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeCrosswalk/NewCode_num_5.txt',header=T)
NewCode = crosswalk$NewCode[match(row_idx, as.numeric(crosswalk$NewCode_num)-100000000)]
UKB_New_500vec = cbind(NewCode,UKB_New_500vec)
UKB_New_500vec$NewCode = sub('Org_','',UKB_New_500vec$NewCode)

####MGI embedding ####
# read New_num for embeddings
load("~/EHRmapping/MGI_embedding/MGI_thres10_pmi/CodeEmbedding/New_0days_500vec.RData")
CodeEmbedVec_MGI = CodeEmbedVec
row_idx = as.numeric(row.names(CodeEmbedVec_MGI))
# read embeddings
MGI_New_500vec = read_csv("~/EHRmapping/MGI_embedding/MGI_thres10_pmi/CodeEmbedding/New_0days_500vec.csv")
MGI_New_500vec = cbind(row_idx,MGI_New_500vec)
# read crosswalk file for mapping
crosswalk = read.table('~/EHRmapping/MGI_embedding/MGI_thres10_pmi/CodeCrosswalk/NewCode_num_10.txt',header=T)
NewCode = crosswalk$NewCode[match(row_idx, as.numeric(crosswalk$NewCode_num)-100000000)]
MGI_New_500vec = cbind(NewCode,MGI_New_500vec)
MGI_New_500vec$NewCode = sub('Org_','',MGI_New_500vec$NewCode)

# normalize embeddings
#### should be normalized based on the chosen vector length.

MGI_New_500vec[,3:(embed_length+2)] = MGI_New_500vec[,3:(embed_length+2)]/apply(MGI_New_500vec[,3:(embed_length+2)],1,norm,type="2")
UKB_New_500vec[,3:(embed_length+2)] = UKB_New_500vec[,3:(embed_length+2)]/apply(UKB_New_500vec[,3:(embed_length+2)],1,norm,type="2")

##### select overlapped codes between the two datasets #####
index_code = intersect(MGI_New_500vec$NewCode,UKB_New_500vec$NewCode) #length(2013)
index_code = index_code[!grepl("_",index_code)]

sum(!grepl("_",UKB_New_500vec$NewCode)) # 9395
sum(!grepl("_",MGI_New_500vec$NewCode)) # 25217

# embedding vector length
lap_MGI = MGI_New_500vec[which(MGI_New_500vec$NewCode %in% index_code), c(1:(2+embed_length))]
lap_UKB = UKB_New_500vec[which(UKB_New_500vec$NewCode %in% index_code), c(1:(2+embed_length))]

# match the rows
lap_MGI = lap_MGI[order(lap_MGI$NewCode),]
lap_UKB = lap_UKB[order(lap_UKB$NewCode),]

# select the embedding vectors
MGI_map = as.matrix(lap_MGI[,(2+1):(2+embed_length)])
UKB_map = as.matrix(lap_UKB[,(2+1):(2+embed_length)])

### datasets for mapping 
UKB_all = UKB_New_500vec
MGI_all = MGI_New_500vec

y = as.matrix(UKB_all[,c(2+(1:embed_length))])
x = as.matrix(MGI_all[,c(2+(1:embed_length))])
rownames(y) = UKB_all$NewCode
rownames(x) = MGI_all$NewCode

###################################
########### Spherical regression ##
###################################
beta_hat = gradient_update_nogrp(X=MGI_map,Y=UKB_map,alpha=1,convergence=1e-4)

coss <- mt_method(x,y,beta_hat)
rownames(coss)=MGI_all$NewCode
colnames(coss)=UKB_all$NewCode

########## For the target group ########
icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]

##### add phecode to embeddings #####
MGI_phecode_round = icdcm_phecode$phecode_round[match(MGI_all$NewCode, icdcm_phecode$ICD_dx)]
UKB_phecode_round = icd_phecode$phecode_round[match(UKB_all$NewCode, icd_phecode$ICD_dx)]
MGI_NewCode_phecode = na.omit(data.frame("NewCode" = MGI_all$NewCode, "phecode_round" = MGI_phecode_round))
UKB_NewCode_phecode = na.omit(data.frame("NewCode" = UKB_all$NewCode, "phecode_round" = UKB_phecode_round))

length(unique(MGI_NewCode_phecode$phecode_round) )
length(unique(UKB_NewCode_phecode$phecode_round) )

overlapped_phecodes = intersect(unique(MGI_NewCode_phecode$phecode_round),
                                unique(UKB_NewCode_phecode$phecode_round)) # 552

results = NULL

overlapped_phecodes = "401"

phecode = "401"
UKB_phecode_ICDx = UKB_NewCode_phecode$NewCode[which(UKB_NewCode_phecode$phecode_round==phecode)]
MGI_phecode_ICDx = MGI_NewCode_phecode$NewCode[which(MGI_NewCode_phecode$phecode_round==phecode)]
# UKB_phecode_ICD9 <- UKB_phecode_ICDx[!grepl("([A-Z][0-9]{2})",UKB_phecode_ICDx)]
UKB_phecode_ICD10 <- UKB_phecode_ICDx[grepl("([A-Z][0-9]{2})",UKB_phecode_ICDx)]
# MGI_phecode_ICD9 <- MGI_phecode_ICDx[grepl("(^[0-9]{2}|[V][0-9]{2})",MGI_phecode_ICDx)]
MGI_phecode_ICD10 <- MGI_phecode_ICDx[grepl("([A-Z][0-9]{2})",MGI_phecode_ICDx)]
MGI_phecode_ICD10 <- MGI_phecode_ICD10[!grepl("(^[V][0-9]{2})",MGI_phecode_ICD10)]
# MGI_phecode_ICD10 <- MGI_phecode_ICD10[!grepl("(^[N][0-9]{2})",MGI_phecode_ICD10)]

##### after harmonization ####
coss_phecode = coss[rownames(coss)%in%MGI_phecode_ICD10,colnames(coss)%in%UKB_phecode_ICD10, drop=FALSE]
# set the frequency
MGI_freq_ICD10_noDuplicates = MGI_freq_all_noDuplicates$freq[match(MGI_phecode_ICD10,MGI_freq_all_noDuplicates$icd10)]
UKB_freq_ICD10_noDuplicates = UKB_freq_all_noDuplicates$freq[match(UKB_phecode_ICD10,UKB_freq_all_noDuplicates$icd10)]
MGI_freq_ICD10_norm_noDuplicates = MGI_freq_ICD10_noDuplicates/sum(MGI_freq_ICD10_noDuplicates)
UKB_freq_ICD10_norm_noDuplicates = UKB_freq_ICD10_noDuplicates/sum(UKB_freq_ICD10_noDuplicates)

MGI_freq_ICD10_withDuplicates = MGI_freq_all_withDuplicates$freq[match(MGI_phecode_ICD10,MGI_freq_all_withDuplicates$icd10)]
UKB_freq_ICD10_withDuplicates = UKB_freq_all_withDuplicates$freq[match(UKB_phecode_ICD10,UKB_freq_all_withDuplicates$icd10)]
MGI_freq_ICD10_norm_withDuplicates = MGI_freq_ICD10_withDuplicates/sum(MGI_freq_ICD10_withDuplicates)
UKB_freq_ICD10_norm_withDuplicates = UKB_freq_ICD10_withDuplicates/sum(UKB_freq_ICD10_withDuplicates)

# using the original cosine matrix
# make the table
# ids
UKB_ids = unique(UKB_ICD10$id)
MGI_ids = unique(MGI_ICD10$id)
# UKB
UKB_data = NULL
for(j in 1:length(UKB_phecode_ICD10)){
  code = UKB_phecode_ICD10[j]
  temp = UKB_ICD10[which(UKB_ICD10$icd10==code),]
  UKB_temp = as.numeric(UKB_ids %in% temp$id)
  UKB_data = cbind(UKB_data, UKB_temp)
}
colnames(UKB_data) = UKB_phecode_ICD10
# MGI
MGI_data = NULL
for(j in 1:length(MGI_phecode_ICD10)){
  code = MGI_phecode_ICD10[j]
  temp = MGI_ICD10[which(MGI_ICD10$icd10==code),]
  MGI_temp = as.numeric(MGI_ids %in% temp$id)
  MGI_data = cbind(MGI_data, MGI_temp)
}
colnames(MGI_data) = MGI_phecode_ICD10

zeroEntry_MGI_all = mean(apply(MGI_data,1,sum)==0)
zerEntry_UKB_all = mean(apply(UKB_data,1,sum)==0)

#### original cosine ###
Gamma_hat = apply(coss_phecode,2,as.numeric)
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.orgCos.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.orgCos = paste0(auc.orgCos.est[2]," (", auc.orgCos.est[1],", ",auc.orgCos.est[3], ") ")

#### top 1 mapping ####
Gamma_hat = t(apply(coss_phecode,1,function(x){ (x==max(x))*x }))
if(ncol(coss_phecode)==1) {Gamma_hat = t(Gamma_hat)}
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.orgCos.top1.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.orgCos.top1 = paste0(auc.orgCos.top1.est[2]," (", auc.orgCos.top1.est[1],", ",auc.orgCos.top1.est[3], ") ")

#### after thresholding #####
thresholds = seq(0.1,0.7,by=0.01)
UKB_embedding_ICD10target = y[which(rownames(y) %in% UKB_phecode_ICD10),]
MGI_embedding_ICD10target = x[which(rownames(x) %in% MGI_phecode_ICD10),]
MGI_ICD10_beta_hat = MGI_embedding_ICD10target %*% beta_hat
loss.evaluate = coss_cv(y_hat = MGI_ICD10_beta_hat,
                        y = UKB_embedding_ICD10target,
                        thres = thresholds,
                        marginal=FALSE, UKB_freq_norm=MGI_freq_ICD10_norm, MGI_freq_norm=UKB_freq_ICD10_norm,
                        nfold=20,
                        leaveOneOut=FALSE)
loss.results = loss.evaluate$loss
(cutoff.ICD10_SR = loss.evaluate$cutoff.min)

coss_org_truncated = coss_phecode
coss_org_truncated[which(coss_org_truncated<cutoff.ICD10_SR)]=0

Gamma_hat = coss_org_truncated
if(ncol(coss_phecode)==1) {Gamma_hat = t(Gamma_hat)}
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.orgCos.threshold.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.orgCos.threshold = paste0(auc.orgCos.threshold.est[2]," (", auc.orgCos.threshold.est[1],", ",auc.orgCos.threshold.est[3], ") ")

##### incorporating the marginal frequency:no Duplicates ######
coss_refined = refine_map_mt(map_mt = coss_phecode, source = MGI_freq_ICD10_norm_noDuplicates, target = UKB_freq_ICD10_norm_noDuplicates)

Gamma_hat = apply(coss_refined,2,as.numeric)
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.refinedCos.noDuplicates.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.refinedCos.noDuplicates = paste0(auc.refinedCos.noDuplicates.est[2]," (", 
                        auc.refinedCos.noDuplicates.est[1],", ",
                        auc.refinedCos.noDuplicates.est[3], ") ")

#### top 1 mapping ####
Gamma_hat = t(apply(coss_phecode,1,function(x){ (x==max(x))*x }))
if(ncol(coss_phecode)==1) {Gamma_hat = t(Gamma_hat)}
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)

# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.refinedCos.noDuplicates.top1.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.refinedCos.noDuplicates.top1 = paste0(auc.refinedCos.noDuplicates.top1.est[2]," (", 
                                          auc.refinedCos.noDuplicates.top1.est[1],", ",
                                          auc.refinedCos.noDuplicates.top1.est[3], ") ")

###### after thresholding #####
thresholds = seq(0.1,0.7,by=0.01)
MGI_ICD10_beta_hat = MGI_embedding_ICD10target %*% beta_hat
loss.evaluate = coss_cv(y_hat = MGI_ICD10_beta_hat,
                        y = UKB_embedding_ICD10target,
                        thres = thresholds,
                        marginal=TRUE, UKB_freq_norm=MGI_freq_ICD10_norm_noDuplicates, 
                        MGI_freq_norm=UKB_freq_ICD10_norm_noDuplicates,
                        nfold=20,
                        leaveOneOut=FALSE)
loss.results = loss.evaluate$loss
(cutoff.ICD10_SR = loss.evaluate$cutoff.min)

coss_refined_truncated = coss_refined
coss_refined_truncated[which(coss_refined_truncated<cutoff.ICD10_SR)]=0

Gamma_hat = coss_refined_truncated
if(ncol(coss_phecode)==1) {Gamma_hat = t(Gamma_hat)}
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.refinedCos.noDuplicates.threshold.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.refinedCos.noDuplicates.threshold = paste0(auc.refinedCos.noDuplicates.threshold.est[2]," (", 
                                               auc.refinedCos.noDuplicates.threshold.est[1],", ",
                                               auc.refinedCos.noDuplicates.threshold.est[3], ") ")


##### incorporating the marginal frequency:with Duplicates ######
coss_refined = refine_map_mt(map_mt = coss_phecode, 
                             source = MGI_freq_ICD10_norm_withDuplicates, 
                             target = UKB_freq_ICD10_norm_withDuplicates)

Gamma_hat = apply(coss_refined,2,as.numeric)
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.refinedCos.withDuplicates.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.refinedCos.withDuplicates = paste0(auc.refinedCos.withDuplicates.est[2]," (", 
                                       auc.refinedCos.withDuplicates.est[1],", ",
                                       auc.refinedCos.withDuplicates.est[3], ") ")

#### top 1 mapping ####
Gamma_hat = t(apply(coss_phecode,1,function(x){ (x==max(x))*x }))
if(ncol(coss_phecode)==1) {Gamma_hat = t(Gamma_hat)}
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)

# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.refinedCos.withDuplicates.top1.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.refinedCos.withDuplicates.top1 = paste0(auc.refinedCos.withDuplicates.top1.est[2]," (", 
                                            auc.refinedCos.withDuplicates.top1.est[1],", ",
                                            auc.refinedCos.withDuplicates.top1.est[3], ") ")

###### after thresholding #####
thresholds = seq(0.1,0.7,by=0.01)
MGI_ICD10_beta_hat = MGI_embedding_ICD10target %*% beta_hat
loss.evaluate = coss_cv(y_hat = MGI_ICD10_beta_hat,
                        y = UKB_embedding_ICD10target,
                        thres = thresholds,
                        marginal=TRUE, UKB_freq_norm=MGI_freq_ICD10_norm_withDuplicates, 
                        MGI_freq_norm=UKB_freq_ICD10_norm_withDuplicates,
                        nfold=20,
                        leaveOneOut=FALSE)
loss.results = loss.evaluate$loss
(cutoff.ICD10_SR = loss.evaluate$cutoff.min)

coss_refined_truncated = coss_refined
coss_refined_truncated[which(coss_refined_truncated<cutoff.ICD10_SR)]=0

Gamma_hat = coss_refined_truncated
if(ncol(coss_phecode)==1) {Gamma_hat = t(Gamma_hat)}
MGI_trans = MGI_data %*% Gamma_hat
data_comb = rbind(cbind(1,MGI_trans),cbind(0,UKB_data))
colnames(data_comb) = c("y",UKB_phecode_ICD10)
data_comb = as.data.frame(data_comb)
# AUC
model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
model = glm(model.formula,data = data_comb,family="binomial")
predicted <- predict(model, data_comb, type="response")
auc.refinedCos.withDuplicates.threshold.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
auc.refinedCos.withDuplicates.threshold = paste0(auc.refinedCos.withDuplicates.threshold.est[2]," (", 
                                                 auc.refinedCos.withDuplicates.threshold.est[1],", ",
                                                 auc.refinedCos.withDuplicates.threshold.est[3], ") ")

##### without harmonization ##############
# number of overlapped codes
codes_overlapped = intersect(MGI_phecode_ICD10, UKB_phecode_ICD10)

if(length(codes_overlapped)==0){
  auc.noHar = 222
}else{
  # UKB
  UKB_data = NULL
  for(j in 1:length(codes_overlapped)){
    code = codes_overlapped[j]
    temp = UKB_ICD10[which(UKB_ICD10$icd10==code),]
    UKB_temp = as.numeric(UKB_ids %in% temp$id)
    UKB_data = cbind(UKB_data, UKB_temp)
  }
  colnames(UKB_data) = codes_overlapped
  # MGI
  MGI_data = NULL
  for(j in 1:length(codes_overlapped)){
    code = codes_overlapped[j]
    temp = MGI_ICD10[which(MGI_ICD10$icd10==code),]
    MGI_temp = as.numeric(MGI_ids %in% temp$id)
    MGI_data = cbind(MGI_data, MGI_temp)
  }
  colnames(MGI_data) = codes_overlapped
  
  data_comb = rbind(cbind(1,MGI_data),cbind(0,UKB_data))
  colnames(data_comb) = c("y",codes_overlapped)
  data_comb = as.data.frame(data_comb)
  # AUC
  model.formula = as.formula(paste0("y~",paste0(codes_overlapped,collapse = "+")))
  model = glm(model.formula,data = data_comb,family="binomial")
  predicted <- predict(model, data_comb, type="response")
  auc.noHar.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
  auc.noHar = paste0(auc.noHar.est[2]," (", auc.noHar.est[1],", ",auc.noHar.est[3], ") ")
  
}

result = c(phecode,
           mean(coss_phecode),
           length(MGI_phecode_ICD10),
           length(UKB_phecode_ICD10),
           length(codes_overlapped),
           auc.orgCos,
           auc.orgCos.top1,
           auc.refinedCos,
           auc.refinedCos.top1,
           auc.noHar,
           zeroEntry_MGI_all,
           zerEntry_UKB_all
)

result0 = data.frame("phecode"=phecode,
                     "auc.woHarmonization" = auc.noHar,
                     "auc.orgCos.top1" = auc.orgCos.top1,
                     "auc.orgCos.threshold" = auc.orgCos.threshold,
                     "auc.refinedCos.noDuplicates.top1" = auc.refinedCos.noDuplicates.top1,
                     "auc.refinedCos.noDuplicates.threshold" = auc.refinedCos.noDuplicates.threshold,
                     "auc.refinedCos.withDuplicates.top1" = auc.refinedCos.withDuplicates.top1,
                     "auc.refinedCos.withDuplicates.threshold" = auc.refinedCos.withDuplicates.threshold
                     )
write.csv(result0,"~/EHRmapping/mapping/results/AUC_target/AUC_ICD10x401_MGItoUKB.csv")





