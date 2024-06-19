library(data.table)
library(dplyr)
library(stringr)
library(SKAT)

source("~/EHRmapping/createUKBphenome/scripts/function.reformatUKB.r")
source("~/EHRmapping/createUKBphenome/scripts/function.harmonizeICD9.r")
source("~/EHRmapping/createUKBphenome/scripts/function.harmonizeICD10.r")

load("~/EHRmapping/mapping/data/demographics.RData")
UKB_demo0 = UKB_demo; MGI_demo0 = MGI_demo

# sourcing the phecode information
icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]

load('~/EHRmapping/comparison/data/UKB_ICD10_noDuplicates.RData')
load('~/EHRmapping/comparison/data/MGI_ICD10_noDuplicates.RData')

# find the block information
colnames(UKB_ICD10)[1:2] = c("id","icd10")
colnames(MGI_ICD10)[1:2] = c("id","icd10")

UKB_ICD10$phecode_round = icd_phecode$phecode_round[match(UKB_ICD10$icd10,icd_phecode$ICD_dx)]
MGI_ICD10$phecode_round = icdcm_phecode$phecode_round[match(MGI_ICD10$icd10,icdcm_phecode$ICD_dx)]

length(unique(UKB_ICD10$icd10)) # 11782
length(unique(MGI_ICD10$icd10)) # 28868

UKB_ICD10_block_unique = unique(na.omit(UKB_ICD10$phecode_round)) # 552
MGI_ICD10_block_unique = unique(na.omit(MGI_ICD10$phecode_round)) # 566
block_overlapped = intersect(UKB_ICD10_block_unique, MGI_ICD10_block_unique) # 549

############ unique blocks in two datasets ########
block_diff = setdiff(UKB_ICD10_block_unique,MGI_ICD10_block_unique)

UKB_unique_block = UKB_ICD10[which(UKB_ICD10$phecode_round %in% block_diff),c(2,3)]
UKB_unique_block = UKB_unique_block[!duplicated(UKB_unique_block),]
counts=NULL
MGI_counts = NULL
MGI_phecodeMap = NULL
for(i in 1:nrow(UKB_unique_block)){
  print(i)
  code = UKB_unique_block$icd10[i]
  UKB_temp = UKB_ICD10[which(UKB_ICD10$icd10==code),]
  count = nrow(UKB_temp)
  counts = c(counts,count)

  # information in MGI
  count_i = nrow(MGI_ICD10[which(MGI_ICD10$icd10==code),c("icd10","phecode_round")])
  MGI_counts = c(MGI_counts,count_i)
  MGI_phecodeMap = c(MGI_phecodeMap, icdcm_phecode$phecode_round[match(code,icdcm_phecode$ICD_dx)])

}
UKB_unique_block = cbind(UKB_unique_block, counts,MGI_counts, MGI_phecodeMap)
colnames(UKB_unique_block) = c("icd10","phecode_icd10","UKB_counts","MGI_counts","phecode_icd10cm")
write.csv(UKB_unique_block,"~/EHRmapping/comparison/results/blockWise_step0_ICD10_uniqueBlock_UKB.csv")

# MGI_unique
(block_diff = setdiff(MGI_ICD10_block_unique,UKB_ICD10_block_unique))
MGI_unique_block = MGI_ICD10[which(MGI_ICD10$phecode_round %in% block_diff),c(2,3)]
MGI_unique_block = MGI_unique_block[!duplicated(MGI_unique_block),]
counts=NULL
UKB_counts = NULL
UKB_phecodeMap = NULL
for(i in 1:nrow(MGI_unique_block)){
  print(i)
  code = MGI_unique_block$icd10[i]
  MGI_temp = MGI_ICD10[which(MGI_ICD10$icd10==code),]
  count = nrow(MGI_temp)
  counts = c(counts,count)

  # information in UKB
  count_i = nrow(UKB_ICD10[which(UKB_ICD10$icd10==code),c("icd10","phecode_round")])
  UKB_counts = c(UKB_counts,count_i)
  UKB_phecodeMap = c(UKB_phecodeMap, icd_phecode$phecode_round[match(code,icd_phecode$ICD_dx)])
}
MGI_unique_block = cbind(MGI_unique_block, counts, UKB_counts, UKB_phecodeMap)
colnames(MGI_unique_block) = c("icd10","phecode_icd10cm","MGI_counts","UKB_counts","phecode_icd10")
write.csv(MGI_unique_block,"~/EHRmapping/comparison/results/blockWise_step0_ICD10_uniqueBlock_MGI.csv")


########## select codes ###########
# should remove nas' in phecodes here.
UKB_ICD10 = na.omit(UKB_ICD10)
MGI_ICD10 = na.omit(MGI_ICD10)
UKB_ICD10_codes = unique(UKB_ICD10$icd10)
MGI_ICD10_codes = unique(MGI_ICD10$icd10)
length(UKB_ICD10_codes) # 6748
length(MGI_ICD10_codes) # 27334

# create the code-block crosswalk file that considers both datasets
UKB_ICD10_overlappedBlock = UKB_ICD10[which(UKB_ICD10$phecode_round %in% block_overlapped),c("icd10","phecode_round")]
MGI_ICD10_overlappedBlock = MGI_ICD10[which(MGI_ICD10$phecode_round %in% block_overlapped),c("icd10","phecode_round")]
UKB_ICD10_overlappedBlock = UKB_ICD10_overlappedBlock[!duplicated(UKB_ICD10_overlappedBlock),]
MGI_ICD10_overlappedBlock = MGI_ICD10_overlappedBlock[!duplicated(MGI_ICD10_overlappedBlock),]

# check the discordant pairs of icd10-phecode mapping.
data_merge = merge(UKB_ICD10_overlappedBlock,MGI_ICD10_overlappedBlock,by="icd10",all=TRUE)
colnames(data_merge) = c("icd10","phecode_r_UKB","phecode_r_MGI")
data_merge = na.omit(data_merge)
discor_indx = which(data_merge$phecode_r_UKB!=data_merge$phecode_r_MGI)
data_discor = data_merge[discor_indx,]
discor_UKB_count=apply(data_discor,1,function(x){sum(UKB_ICD10$icd10==x[1])})
discor_MGI_count=apply(data_discor,1,function(x){sum(MGI_ICD10$icd10==x[1])})
data_discor$UKB_count = discor_UKB_count
data_discor$MGI_count = discor_MGI_count

# write.csv(data_discor,"~/EHRmapping/comparison/results/blockWise_step0_ICD10_discordantPairs.csv")
