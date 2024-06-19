library(data.table)
library(dplyr)
library(stringr)
library(SKAT)

load("~/EHRmapping/mapping/data/demographics.RData")
UKB_demo0 = UKB_demo; MGI_demo0 = MGI_demo

# sourcing the phecode information
icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]

load('~/EHRmapping/comparison/data/UKB_ICD10_noDuplicates.RData')
load('~/EHRmapping/comparison/data/MGI_ICD10_noDuplicates.RData')

# find the block information
colnames(UKB_ICD10) = c("id","icd10")
colnames(MGI_ICD10) = c("id","icd10")

UKB_ICD10$phecode_round = icd_phecode$phecode_round[match(UKB_ICD10$icd10,icd_phecode$ICD_dx)]
MGI_ICD10$phecode_round = icdcm_phecode$phecode_round[match(MGI_ICD10$icd10,icdcm_phecode$ICD_dx)]

length(unique(UKB_ICD10$icd10)) # 11782
length(unique(MGI_ICD10$icd10)) # 28868

UKB_ICD10_block_unique = unique(na.omit(UKB_ICD10$phecode_round)) # 552
MGI_ICD10_block_unique = unique(na.omit(MGI_ICD10$phecode_round)) # 566
block_overlapped = intersect(UKB_ICD10_block_unique, MGI_ICD10_block_unique) # 549

########## select codes ###########
# should remove nas' in phecodes here.
UKB_ICD10 = na.omit(UKB_ICD10)
MGI_ICD10 = na.omit(MGI_ICD10)
UKB_ICD10_codes = unique(UKB_ICD10$icd10)
MGI_ICD10_codes = unique(MGI_ICD10$icd10)
length(UKB_ICD10_codes) 
length(MGI_ICD10_codes) 

# create the code-block crosswalk file that considers both datasets
UKB_ICD10_overlappedBlock = UKB_ICD10[which(UKB_ICD10$phecode_round %in% block_overlapped),c("icd10","phecode_round")]
MGI_ICD10_overlappedBlock = MGI_ICD10[which(MGI_ICD10$phecode_round %in% block_overlapped),c("icd10","phecode_round")]
UKB_ICD10_overlappedBlock = UKB_ICD10_overlappedBlock[!duplicated(UKB_ICD10_overlappedBlock),]
MGI_ICD10_overlappedBlock = MGI_ICD10_overlappedBlock[!duplicated(MGI_ICD10_overlappedBlock),]

# remove discordant codes
data_discor = read.csv("~/EHRmapping/comparison/results/blockWise_step0_ICD10_discordantPairs.csv")

block_code_crosswalk = rbind(UKB_ICD10_overlappedBlock,MGI_ICD10_overlappedBlock)
block_code_crosswalk = block_code_crosswalk[!duplicated(block_code_crosswalk),]
block_code_crosswalk = block_code_crosswalk %>%
  arrange(phecode_round,icd10)
block_code_crosswalk = block_code_crosswalk[!(block_code_crosswalk$icd10 %in% data_discor$icd10),]

dim(block_code_crosswalk) # 30379     2
colnames(block_code_crosswalk) = c("code","block")

#####################
UKB_ICD10_selectedCodes = unique(UKB_ICD10$icd10[which(UKB_ICD10$phecode_round %in% block_overlapped)]) #6739
MGI_ICD10_selectedCodes = unique(MGI_ICD10$icd10[which(MGI_ICD10$phecode_round %in% block_overlapped)]) #26999
codes_joint = union(UKB_ICD10_selectedCodes, MGI_ICD10_selectedCodes) #30310
codes_joint = setdiff(codes_joint,data_discor$icd10) # removed discordant codes from comparison, 30241

# codes_joint = unique(block_code_crosswalk$code)

UKB_ICD10_joint = UKB_ICD10[which(UKB_ICD10$icd10 %in% codes_joint),]
MGI_ICD10_joint = MGI_ICD10[which(MGI_ICD10$icd10 %in% codes_joint),]


# rm(list = setdiff(ls(),c("UKB_ICD10_joint","MGI_ICD10_joint","block_code_crosswalk")))
save.image("~/EHRmapping/comparison/results/blockWise_step1_3_ICD10_prepareData.RData")

