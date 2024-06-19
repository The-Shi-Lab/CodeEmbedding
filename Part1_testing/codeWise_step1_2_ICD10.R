# ###############################
# #### adjusted code-wise comparison #####
# ###############################
# 
library(data.table)
library(dplyr)
library(stringr)
library(dplyr)


load("~/EHRmapping/mapping/data/demographics.RData")
UKB_demo0 = UKB_demo; MGI_demo0 = MGI_demo

# sourcing the phecode information
icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]

############# source data ###############
load('~/EHRmapping/comparison/data/UKB_ICD10_noDuplicates.RData')
load('~/EHRmapping/comparison/data/MGI_ICD10_noDuplicates.RData')

#############################
# original number of codes
length(unique(UKB_ICD10$ICD)) # 3500
length(unique(MGI_ICD10$ICD)) # 11934

# choose the overlapped codes
codes_overlapped = intersect(unique(UKB_ICD10$ICD), unique(MGI_ICD10$ICD)) # 4185

UKB_ICD10_overlapped = UKB_ICD10[which(UKB_ICD10$ICD %in% codes_overlapped),] #[1] 2080965       2
MGI_ICD10_overlapped = MGI_ICD10[which(MGI_ICD10$ICD %in% codes_overlapped),] # [1] 1922451       2
colnames(UKB_ICD10_overlapped) = c("id","ICD")
colnames(MGI_ICD10_overlapped) = c("id","ICD")

# create the block information
UKB_ICD10_overlapped$phecode_round = icd_phecode$phecode_round[match(UKB_ICD10_overlapped$ICD,icd_phecode$ICD_dx)] # using phecode = 1084, using phecode_round = 433
MGI_ICD10_overlapped$phecode_round = icdcm_phecode$phecode_round[match(MGI_ICD10_overlapped$ICD,icdcm_phecode$ICD_dx)]

length(unique(UKB_ICD10_overlapped$ICD))
length(unique(MGI_ICD10_overlapped$ICD))


# create wide format
UKB_id = unique(UKB_ICD10_overlapped$id)
UKB_demo = UKB_demo0[match(UKB_id,UKB_demo0$id),]

UKB_ICD10_wide = matrix(NA, nrow = length(UKB_id), ncol=length(codes_overlapped))
for(j in 1:length(codes_overlapped)){
  print(j)
  code = codes_overlapped[j]
  selectedID = unique(UKB_ICD10_overlapped$id[which(UKB_ICD10_overlapped$ICD == code)])
  UKB_ICD10_wide[,j] = (UKB_id %in% selectedID)*1
}

MGI_id = unique(MGI_ICD10_overlapped$id)
MGI_demo = MGI_demo0[match(MGI_id,MGI_demo0$id),]

MGI_ICD10_wide = matrix(NA, nrow = length(MGI_id), ncol=length(codes_overlapped))
for(j in 1:length(codes_overlapped)){
  print(j)
  code = codes_overlapped[j]
  selectedID = unique(MGI_ICD10_overlapped$id[which(MGI_ICD10_overlapped$ICD == code)])
  MGI_ICD10_wide[,j] = (MGI_id %in% selectedID)*1
}

colnames(UKB_ICD10_wide) = codes_overlapped
colnames(MGI_ICD10_wide) = codes_overlapped


rm(list = setdiff(ls(),c("UKB_demo","MGI_demo","UKB_ICD10_wide","MGI_ICD10_wide")))

save.image("~/EHRmapping/comparison/results/codeWise_step1_1_data_ICD10_step1.RData")

####################################
#### Code-wise comparison 2 ########
####################################
load("~/EHRmapping/comparison/results/codeWise_step1_1_data_ICD10_step1.RData")
all_codes = colnames(UKB_ICD10_wide)
print("start code-wise comparison")

pvalues = rep(NA, length(all_codes))

for(c in 1:length(all_codes)){
  print(c)
  code = all_codes[c]
  UKB_data_code = UKB_ICD10_wide[,which(colnames(UKB_ICD10_wide)==code)]
  MGI_data_code = MGI_ICD10_wide[,which(colnames(MGI_ICD10_wide)==code)]
  UKB_data = cbind(UKB_data_code,UKB_demo[,-1]); colnames(UKB_data)[1] = "code"
  MGI_data = cbind(MGI_data_code,MGI_demo[,-1]); colnames(MGI_data)[1] = "code"
  UKB_data = na.omit(UKB_data)
  MGI_data = na.omit(MGI_data)
  
  y = c(rep(1,nrow(UKB_data)),rep(0,nrow(MGI_data)))
  data_comb = as.data.frame(cbind(y,rbind(UKB_data,MGI_data)))
  glm.model = glm(y~code+age+sex,data=data_comb,family="binomial")
  pvalues[c] = summary(glm.model)$coef[2,4]
  # t.test.model = t.test(x=UKB_freq, y=MGI_freq)
  # pvalues[c] = t.test.model$p.value
}

pvalues.codeWise.df = data.frame("code" = all_codes, "pvalue"=pvalues)

# rm(list = setdiff(ls(),c("pvalues.codeWise.df")))
save.image("~/EHRmapping/comparison/results/codeWise_step1_1_data_ICD10_step2.RData")

# pvalues.codeWise.df


