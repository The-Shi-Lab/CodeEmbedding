library(data.table)
library(dplyr)
library(stringr)
library(SKAT)


load("~/EHRmapping/mapping/data/demographics_ICD9.RData")
UKB_demo0 = UKB_demo; MGI_demo0 = MGI_demo

# sourcing the phecode information
icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]

############# source data ###############
load('~/EHRmapping/comparison/data/UKB_ICD9_noDuplicates.RData')
load('~/EHRmapping/comparison/data/MGI_ICD9_noDuplicates.RData')

# find the block information
colnames(UKB_ICD9) = c("id","icd9")
colnames(MGI_ICD9) = c("id","icd9")

UKB_ICD9$phecode_round = icd_phecode$phecode_round[match(UKB_ICD9$icd9,icd_phecode$ICD_dx)]
MGI_ICD9$phecode_round = icdcm_phecode$phecode_round[match(MGI_ICD9$icd9,icdcm_phecode$ICD_dx)]

length(unique(UKB_ICD9$icd9)) 
length(unique(MGI_ICD9$icd9)) 

UKB_ICD9_block_unique = unique(na.omit(UKB_ICD9$phecode_round)) # 552
MGI_ICD9_block_unique = unique(na.omit(MGI_ICD9$phecode_round)) # 566
block_overlapped = intersect(UKB_ICD9_block_unique, MGI_ICD9_block_unique) # 549

########## select codes ###########
# should remove nas' in phecodes here.
UKB_ICD9 = na.omit(UKB_ICD9)
MGI_ICD9 = na.omit(MGI_ICD9)
UKB_ICD9_codes = unique(UKB_ICD9$icd9)
MGI_ICD9_codes = unique(MGI_ICD9$icd9)
length(UKB_ICD9_codes) 
length(MGI_ICD9_codes) 

# create the code-block crosswalk file that considers both datasets
UKB_ICD9_overlappedBlock = UKB_ICD9[which(UKB_ICD9$phecode_round %in% block_overlapped),c("icd9","phecode_round")]
MGI_ICD9_overlappedBlock = MGI_ICD9[which(MGI_ICD9$phecode_round %in% block_overlapped),c("icd9","phecode_round")]
UKB_ICD9_overlappedBlock = UKB_ICD9_overlappedBlock[!duplicated(UKB_ICD9_overlappedBlock),]
MGI_ICD9_overlappedBlock = MGI_ICD9_overlappedBlock[!duplicated(MGI_ICD9_overlappedBlock),]

block_code_crosswalk = rbind(UKB_ICD9_overlappedBlock,MGI_ICD9_overlappedBlock)
block_code_crosswalk = block_code_crosswalk[!duplicated(block_code_crosswalk),]
block_code_crosswalk = block_code_crosswalk %>%
  arrange(phecode_round,icd9)
dim(block_code_crosswalk) 
colnames(block_code_crosswalk) = c("code","block")

#####################
UKB_ICD9_selectedCodes = unique(UKB_ICD9$icd9[which(UKB_ICD9$phecode_round %in% block_overlapped)]) #6739
MGI_ICD9_selectedCodes = unique(MGI_ICD9$icd9[which(MGI_ICD9$phecode_round %in% block_overlapped)]) #26999
codes_joint = union(UKB_ICD9_selectedCodes, MGI_ICD9_selectedCodes) #30310

UKB_ICD9_joint = UKB_ICD9[which(UKB_ICD9$icd9 %in% codes_joint),]
MGI_ICD9_joint = MGI_ICD9[which(MGI_ICD9$icd9 %in% codes_joint),]


# rm(list = setdiff(ls(),c("UKB_ICD9_joint","MGI_ICD9_joint","block_code_crosswalk")))
save.image("~/EHRmapping/comparison/results/blockWise_step1_1_ICD9_prepareData.RData")
