library(tidyverse)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(pROC)

# read embeddings
UKB_optimalWindow = 129
UKB_optimalVecLen = 500
MGI_optimalWindow = 0
MGI_optimalVecLen = 500
embed_length = 250
MGI_embed_duplicates = TRUE
MGI_embedding_file = "MGI_thres10_mainDiagnosis_inPatient"
# MGI_embed_duplicates = FALSE
# MGI_embedding_file = "MGI_thres10_mainDiagnosis_noDuplicates_inPatient"
phecode = "401" # select target group


outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"

source("~/EHRmapping/mapping2/sourceFun/function_mapping.R")
source("~/EHRmapping/mapping2/sourceFun/refine_mapping.R")

####### Step 1.1: read embeddings #######
#### UKB embedding ####
UKB_row_idx_loc = paste0("/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeEmbedding/New_",UKB_optimalWindow,"days_",UKB_optimalVecLen,"vec.RData")
load(UKB_row_idx_loc)
UKB_row_idx = as.numeric(row.names(CodeEmbedVec))
UKB_embeddings = read_csv(paste0("/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeEmbedding/New_",UKB_optimalWindow,"days_",UKB_optimalVecLen,"vec.csv"))
UKB_crosswalk = read.table('/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeCrosswalk/NewCode_num_5.txt',header=T)
NewCode = UKB_crosswalk$NewCode[match(UKB_row_idx, as.numeric(UKB_crosswalk$NewCode_num)-100000000)]
UKB_embeddings = cbind(NewCode,UKB_embeddings)
UKB_embeddings$NewCode = sub('Org_','',UKB_embeddings$NewCode)
head(UKB_embeddings[,1:5])

#### MGI embedding ####
MGI_row_idx_loc = paste0("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/",MGI_embedding_file,"/CodeEmbedding/New_",MGI_optimalWindow,"days_",MGI_optimalVecLen,"vec.RData")
load(MGI_row_idx_loc)
MGI_row_idx = as.numeric(row.names(CodeEmbedVec))
MGI_embeddings = read_csv(paste0("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/",MGI_embedding_file,"/CodeEmbedding/New_",MGI_optimalWindow,"days_",MGI_optimalVecLen,"vec.csv"))
MGI_crosswalk = read.table(paste0("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/",MGI_embedding_file,"/CodeCrosswalk/NewCode_num_10.txt"),header=T)
NewCode = MGI_crosswalk$NewCode[match(MGI_row_idx, as.numeric(MGI_crosswalk$NewCode_num)-100000000)]
MGI_embeddings = cbind(NewCode,MGI_embeddings)
MGI_embeddings$NewCode = sub('Org_','',MGI_embeddings$NewCode)
head(MGI_embeddings[,1:5])

###### Step 1.2: normalize embeddings ######
UKB_embeddings[,2:(embed_length+1)] = UKB_embeddings[,2:(embed_length+1)]/apply(UKB_embeddings[,2:(embed_length+1)],1,norm,type="2")
MGI_embeddings[,2:(embed_length+1)] = MGI_embeddings[,2:(embed_length+1)]/apply(MGI_embeddings[,2:(embed_length+1)],1,norm,type="2")


###### Step 1.3: space alignment ######
index_code = intersect(UKB_embeddings$NewCode,MGI_embeddings$NewCode) #length(2013)
# index_code = index_code[!grepl("_",index_code)] # should we also align codes that's been grouped?

sum(!grepl("_",UKB_embeddings$NewCode)) # 9395
sum(!grepl("_",MGI_embeddings$NewCode)) # 25217

# embedding vector length
lap_MGI = MGI_embeddings[which(MGI_embeddings$NewCode %in% index_code), c(1:(1+embed_length))]
lap_UKB = UKB_embeddings[which(UKB_embeddings$NewCode %in% index_code), c(1:(1+embed_length))]

# match the rows
lap_MGI = as.matrix(lap_MGI[order(lap_MGI$NewCode),-1])
lap_UKB = as.matrix(lap_UKB[order(lap_UKB$NewCode),-1])

# estimate the translation matrix
beta_hat = gradient_update_nogrp(X=lap_MGI,Y=lap_UKB,alpha=1,convergence=1e-4)

####### Step 1.4: calculate the cosine value #######
UKB_all = as.matrix(UKB_embeddings[,2:(embed_length+1)])
MGI_all = as.matrix(MGI_embeddings[,2:(embed_length+1)])
rownames(UKB_all) = UKB_embeddings$NewCode
rownames(MGI_all) = MGI_embeddings$NewCode

icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]
UKB_selectedICD = find_ICD(phecode=phecode, mappingFile=icd_phecode,poolICD = UKB_embeddings$NewCode, ICDversion="ICD10")
MGI_selectedICD = find_ICD(phecode=phecode, mappingFile=icdcm_phecode,poolICD = MGI_embeddings$NewCode, ICDversion = "ICD10")

UKB_embeddings_map = UKB_all[match(UKB_selectedICD, rownames(UKB_all)),]
MGI_embeddings_map = MGI_all[match(MGI_selectedICD, rownames(MGI_all)),]

(coss = mt_method(MGI_embeddings_map, UKB_embeddings_map,beta_hat))


if(MGI_embed_duplicates){
    outfile = paste0(outpath,"phecode_",phecode,"_orgCosine_withDuplicates.csv")
    outfile_Rdata = paste0(outpath,"phecode_",phecode,"_orgCosine_withDuplicates.RData")
} else {
    outfile = paste0(outpath,"phecode_",phecode,"_orgCosine_noDuplicates.csv")
    outfile_Rdata = paste0(outpath,"phecode_",phecode,"_orgCosine_noDuplicates.RData")}


write.csv(coss, outfile, row.names = TRUE)
# save UKB_embeddings_map and MGI_embeddings_map to outfile_Rdata
save(beta_hat, UKB_embeddings_map, MGI_embeddings_map, file = outfile_Rdata)
