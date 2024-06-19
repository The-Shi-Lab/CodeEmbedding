library(data.table)
library(dplyr)

############ UKB data ##########
load("~/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/Step1_beforeRareCode.RData")
UKB_ICD10_1 = all_icd10_time_copy[,c("id","code")]
UKB_ICD10 = UKB_ICD10_1[!duplicated(UKB_ICD10_1),];colnames(UKB_ICD10) = c("id","ICD")
dim(UKB_ICD10)

UKB_ICD9_1 = all_icd9_time_copy[,c("id","code")]
UKB_ICD9 = UKB_ICD9_1[!duplicated(UKB_ICD9_1),];colnames(UKB_ICD9) = c("id","ICD")
dim(UKB_ICD9)

save(UKB_ICD10, file = "~/EHRmapping/comparison/data/UKB_ICD10_noDuplicates.RData")
save(UKB_ICD9, file = "~/EHRmapping/comparison/data/UKB_ICD9_noDuplicates.RData")


### MGI ####
load("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis_noDuplicates/Step1_beforeRareCodes.RData")
dim(dx)
dx = dx[,c("Deid_ID","DxCode","Lexicon")]
dx_noDuplicates = dx[!duplicated(dx),]
dim(dx_noDuplicates)
MGI_ICD10 = dx_noDuplicates[which(dx_noDuplicates$Lexicon == "ICD10"),c("Deid_ID","DxCode")];colnames(MGI_ICD10) = c("id","ICD")
MGI_ICD9 = dx_noDuplicates[which(dx_noDuplicates$Lexicon == "ICD9"),c("Deid_ID","DxCode")];colnames(MGI_ICD9) = c("id","ICD")
save(MGI_ICD10, file = "~/EHRmapping/comparison/data/MGI_ICD10_noDuplicates.RData")
save(MGI_ICD9, file = "~/EHRmapping/comparison/data/MGI_ICD9_noDuplicates.RData")
