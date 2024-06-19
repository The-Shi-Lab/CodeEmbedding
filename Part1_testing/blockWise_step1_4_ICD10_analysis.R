library(data.table)
library(dplyr)
library(stringr)
library(SKAT)


########## block-wise comparison #########
load("~/EHRmapping/comparison/results/blockWise_step1_3_ICD10_prepareData.RData")

all_codes = block_code_crosswalk$code
all_blocks = unique(block_code_crosswalk$block)

print("start burden test")

pvalues = rep(NA, length(all_blocks))
UKB_id = unique(UKB_ICD10_joint$id)
MGI_id = unique(MGI_ICD10_joint$id)
UKB_demo = UKB_demo0[match(UKB_id,UKB_demo0$id),]
MGI_demo = MGI_demo0[match(MGI_id,MGI_demo0$id),]


for(b in 1:length(all_blocks)){
  print(b)
  block = all_blocks[b]
  codes = block_code_crosswalk$code[which(block_code_crosswalk$block==block)]
  
  UKB_ICD10_wide = MGI_ICD10_wide = NULL
  for(j in 1:length(codes)){
    code = codes[j]
    
    selectedID = unique(MGI_ICD10_joint$id[which(MGI_ICD10_joint$icd10 == code)])
    temp = (MGI_id %in% selectedID)*1
    MGI_ICD10_wide = cbind(MGI_ICD10_wide, temp)
    
    selectedID = unique(UKB_ICD10_joint$id[which(UKB_ICD10_joint$icd10 == code)])
    temp = (UKB_id %in% selectedID)*1
    UKB_ICD10_wide = cbind(UKB_ICD10_wide, temp)
    
  }
  UKB_data = as.matrix(UKB_ICD10_wide)
  MGI_data = as.matrix(MGI_ICD10_wide)
  UKB_freq = apply(UKB_data,1,sum)
  MGI_freq = apply(MGI_data,1,sum)
  
  code = c(UKB_freq,MGI_freq)
  y = c(rep(1,nrow(UKB_data)), rep(0, nrow(MGI_data)))
  data_comb = cbind(y, cbind(code,rbind(UKB_demo, MGI_demo)))
  
  glm.model = glm(y~code+age+sex, data=data_comb, family="binomial")
  pvalues[b] = summary(glm.model)$coef[2,4]
  # t.test.model = t.test(UKB_freq, MGI_freq)
  # pvalues[b] = t.test.model$p.value
}

pvalues.burden.df = data.frame("block" = all_blocks,"pvalue"=pvalues)

######## SKAT test #####

pvalues = rep(NA, length(all_blocks))
demo = rbind(UKB_demo, MGI_demo)
selection = c(rep(1,length(UKB_id)), rep(0,length(MGI_id)))
df_temp = cbind(selection,demo)
# df_temp = na.omit(cbind(selection,demo))
obj<-SKAT_Null_Model(selection~age+sex+race+comor, data= df_temp, out_type="D",Adjustment=F)


for(b in 1:length(all_blocks)){
  print(b)
  block = all_blocks[b]
  codes = block_code_crosswalk$code[which(block_code_crosswalk$block==block)]
  UKB_ICD10_wide = MGI_ICD10_wide = NULL
  for(j in 1:length(codes)){
    code = codes[j]
    
    selectedID = unique(MGI_ICD10_joint$id[which(MGI_ICD10_joint$icd10 == code)])
    temp = (MGI_id %in% selectedID)*1
    MGI_ICD10_wide = cbind(MGI_ICD10_wide, temp)
    
    selectedID = unique(UKB_ICD10_joint$id[which(UKB_ICD10_joint$icd10 == code)])
    temp = (UKB_id %in% selectedID)*1
    UKB_ICD10_wide = cbind(UKB_ICD10_wide, temp)
    
  }
  UKB_data = as.matrix(UKB_ICD10_wide)
  MGI_data = as.matrix(MGI_ICD10_wide)
  tmp.visit = as.matrix(rbind(UKB_data,MGI_data))
  weight.visit = 1/colMeans(tmp.visit,na.rm=TRUE)
  weight.visit = weight.visit/sum(weight.visit)
  model= SKAT(tmp.visit, obj, weights= weight.visit, is_check_genotype = F)
  pvalues[b] = model$p.value
}

pvalues.skat.df = data.frame("block" = all_blocks,"pvalue"=pvalues)

rm(list = setdiff(ls(),c("pvalues.burden.df","pvalues.skat.df","block_code_crosswalk")))
save.image("~/EHRmapping/comparison/results/blockWise_step1_4_ICD10_analysiss.RData")

