
library(ggplot2)
library(reshape2)
##### import pvalues for code-wise comparison #####
load("~/EHRmapping/comparison/results/codeWise_step1_1_data_ICD10_step2.RData")
rm(list=setdiff(ls(),"pvalues.codeWise.df"))

#### import pvalues for block-wise comparison ####
load("~/EHRmapping/comparison/results/blockWise_step1_4_ICD10_analysiss.RData")

pvalues.burden = pvalues.burden.df
pvalues.skat = pvalues.skat.df

##### 
code_block = block_code_crosswalk
all_blocks = unique(code_block$block)
all_codes = unique(block_code_crosswalk$code)

block_new = data.frame("block" = all_blocks,"block_new" = 1:length(all_blocks))
code_block$block_new = block_new$block_new[match(code_block$block, block_new$block)]
pvalues_results = data.frame("codes" = all_codes, 
                             "pvalues_codes"=pvalues.codeWise.df$pvalue[match(all_codes,pvalues.codeWise.df$code)],
                             "block"= code_block$block[match(all_codes,code_block$code)],
                             "block_new"= code_block$block_new[match(all_codes,code_block$code)])
pvalues_results2 = pvalues_results[order(pvalues_results$block_new),]
pvalues_results2 = na.omit(pvalues_results2)

pvalues.burden.df = data.frame('block'=all_blocks,pvalues=pvalues.burden$pvalue[match(all_blocks,pvalues.burden$block)])
pvalues.burden.df$block_new = block_new$block_new[match(pvalues.burden.df$block, block_new$block)]
pvalues.skat.df = data.frame('block'=all_blocks,pvalues=pvalues.skat$pvalue[match(all_blocks,pvalues.skat$block)])
pvalues.skat.df$block_new = block_new$block_new[match(pvalues.skat.df$block, block_new$block)]

pvalues.block = pvalues.burden.df
pvalues.block$SKAT  = pvalues.skat.df$pvalues[match(pvalues.block$block,pvalues.skat.df$block)]
colnames(pvalues.block)  = c("block","Burden","block_new","SKAT")


##### make the table #####
load("~/EHRmapping/comparison/results/blockWise_step1_3_ICD10_prepareData.RData")

# function for unique blocks in data2
uniqueSigBlock_fun = function(p.threshold, data1, data2){
    ## block-burden
    sig.data1.ind = data1$block[which(data1$pvalue<=p.threshold)]
    length(sig.data1.ind) 

    ## block-skat
    sig.data2.ind = data2$block[which(data2$pvalue<=p.threshold)]
    length(sig.data2.ind) 

    block_indx = setdiff(sig.data2.ind, sig.data1.ind)
    length(block_indx)
    
    print(paste0("The number of unique block is ",length(block_indx)))
    if(length(block_indx)>20){
        break
    }

    df=NULL
    UKB_id = unique(UKB_ICD10_joint$id)
    MGI_id = unique(MGI_ICD10_joint$id)
    for(i in 1:length(block_indx)){
    print(i)
    block=block_indx[i]
    p.data1 = data1$pvalues[which(data1$block==block)]
    p.data2 = data2$pvalues[which(data2$block==block)]
    codes = block_code_crosswalk$code[which(block_code_crosswalk$block==block)]
    
    MGI_ICD10_wide = UKB_ICD10_wide = NULL
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
    
    df_temp = data.frame("UKB"=apply(UKB_data,2,sum), "MGI"= apply(MGI_data,2,sum))
    rownames(df_temp) = codes
    total = data.frame("subtotal"=apply(df_temp,2,sum))
    df_temp = rbind(df_temp,t(total))
    UKB_norm = df_temp$UKB/nrow(UKB_data)
    MGI_norm = df_temp$MGI/nrow(MGI_data)
    df_temp = cbind(df_temp, UKB_norm, MGI_norm)
    df_temp = cbind(block,p.data2,p.data1,df_temp)
    
    df= rbind(df,df_temp)
    }
    return(df)
}

df = uniqueSigBlock_fun(p.threshold=0.05/length(all_blocks), data1 = pvalues.burden.df, data2 = pvalues.skat.df)
write.csv(df,"~/EHRmapping/comparison/results/step3_2_1_ICD10_skat_vs_burden_pBonferroni.csv")

df = uniqueSigBlock_fun(p.threshold=0, data1 = pvalues.burden.df, data2 = pvalues.skat.df)
write.csv(df,"~/EHRmapping/comparison/results/step3_2_1_ICD10_skat_vs_burden_p0.csv")

# df = uniqueSigBlock_fun(p.threshold=0.05/length(all_blocks), data1 = pvalues.skat.df, data2 = pvalues.burden.df)
# write.csv(df,"~/EHRmapping/comparison/results/step3_2_1_ICD10_burden_vs_skat_pBonferroni.csv")

# df = uniqueSigBlock_fun(p.threshold=0, data1 = pvalues.burden.df, data2 = pvalues.skat.df)
# write.csv(df,"~/EHRmapping/comparison/results/step3_2_1_ICD10_burden_vs_skat_p0.csv")

### unique in SKAT, Bonferroni
data = read.csv("~/EHRmapping/comparison/results/step3_2_1_ICD10_skat_vs_burden_pBonferroni.csv")
data1 = data[grepl("subtotal",data$X),]

# UKB& MGI explanation
UKB_phecode = read.csv('/net/snowwhite/home/jiacong/EHRmapping/comparison/data/phecode_icd10.csv')
UKB_phecode = na.omit(UKB_phecode)
colnames(UKB_phecode)[1:3] = c("icd10","icd10.string","phecode")
MGI_phecode = read.csv('/net/snowwhite/home/jiacong/EHRmapping/comparison/data/Phecode_map_v1_2_icd10cm_beta.csv')
colnames(MGI_phecode)[1:3] = c("icd10","icd10.string","phecode")

data2 = data1 %>%
  mutate(phenotype_WHO = UKB_phecode$Phenotype[match(data1$block,UKB_phecode$phecode)]) %>%
  mutate(phenotype_US = MGI_phecode$phecode_str[match(data1$block, MGI_phecode$phecode)]) %>%
  select("block","phenotype_WHO","phenotype_US","p.data2","p.data1","UKB","MGI","UKB_norm","MGI_norm")
colnames(data2) = c("Phecode","Phenotype_WHO","Phenotype_US","p_skat","p_burden","UKB_count","MGI_count","UKB_normedCount","MGI_normedCount")
data2
write.csv(data2, "~/EHRmapping/comparison/results/step3_2_1_ICD10_skat_vs_burden_pBonferroni_detail.csv")

### unique in SKAT, p=0
data = read.csv("~/EHRmapping/comparison/results/step3_2_1_ICD10_skat_vs_burden_p0.csv")
data1 = data[grepl("subtotal",data$X),]

# UKB& MGI explanation
UKB_phecode = read.csv('/net/snowwhite/home/jiacong/EHRmapping/comparison/data/phecode_icd10.csv')
UKB_phecode = na.omit(UKB_phecode)
colnames(UKB_phecode)[1:3] = c("icd10","icd10.string","phecode")
MGI_phecode = read.csv('/net/snowwhite/home/jiacong/EHRmapping/comparison/data/Phecode_map_v1_2_icd10cm_beta.csv')
colnames(MGI_phecode)[1:3] = c("icd10","icd10.string","phecode")

data2 = data1 %>%
  mutate(phenotype_WHO = UKB_phecode$Phenotype[match(data1$block,UKB_phecode$phecode)]) %>%
  mutate(phenotype_US = MGI_phecode$phecode_str[match(data1$block, MGI_phecode$phecode)]) %>%
  select("block","phenotype_WHO","phenotype_US","p.data2","p.data1","UKB","MGI","UKB_norm","MGI_norm")
colnames(data2) = c("Phecode","Phenotype_WHO","Phenotype_US","p_skat","p_burden","UKB_count","MGI_count","UKB_normedCount","MGI_normedCount")
data2
write.csv(data2, "~/EHRmapping/comparison/results/step3_2_1_ICD10_skat_vs_burden_p0_detail.csv")


# ### unique in Burden
# data = read.csv("~/EHRmapping/mapping/results/blockCompare/ICD10_blockCompareAdjusted_uniqueBurden_raw.csv")
# data1 = data[grepl("subtotal",data$X),]

# data2 = data1 %>%
#   mutate(phenotype_WHO = UKB_phecode$Phenotype[match(data1$block,UKB_phecode$phecode)]) %>%
#   mutate(phenotype_US = MGI_phecode$phecode_str[match(data1$block, MGI_phecode$phecode)]) %>%
#   select("block","phenotype_WHO","phenotype_US","p.skat","p.burden","UKB","MGI","UKB_norm","MGI_norm")
# colnames(data2) = c("Phecode","Phenotype_WHO","Phenotype_US","p_skat","p_burden","UKB_count","MGI_count","UKB_normedCount","MGI_normedCount")
# write.csv(data2, "~/EHRmapping/comparison/results/step3_2_1_ICD10block_pBon_uniqueBurden.csv")


