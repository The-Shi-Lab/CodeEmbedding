library(data.table)
library(tidyverse)
library(pROC)
source("~/EHRmapping/mapping2/sourceFun/function_mapping.R")
load("~/EHRmapping/mapping2/data/designMatrix_org.RData")
dim(UKB_ICD)
dim(MGI_ICD)


########################
# read the cosine matrix
MGI_embed_duplicates = FALSE
outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"
cosineType = "orgCosine" 
# cosineType = "refinedCosine"
# sparsify method 
sparsify = "top1"
# sparsify = "noSparsify"
# sparsify = "dataDrivenThreshold"

########## read the embedding file ##############
UKB_optimalWindow = 129
UKB_optimalVecLen = 500
MGI_optimalWindow = 0
MGI_optimalVecLen = 500
embed_length = 500
# embed_length = 250
MGI_embed_duplicates = TRUE
MGI_embedding_file = "MGI_thres10_mainDiagnosis"
# MGI_embed_duplicates = FALSE
# MGI_embedding_file = "MGI_thres10_mainDiagnosis_noDuplicates"
phecode = "401" # select target group

UKB_row_idx_loc = paste0("/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeEmbedding/New_",UKB_optimalWindow,"days_",UKB_optimalVecLen,"vec.RData")
load(UKB_row_idx_loc)
UKB_row_idx = as.numeric(row.names(CodeEmbedVec))
UKB_embeddings = read_csv(paste0("/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeEmbedding/New_",UKB_optimalWindow,"days_",UKB_optimalVecLen,"vec.csv"))
UKB_crosswalk = read.table('/net/snowwhite/home/jiacong/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeCrosswalk/NewCode_num_5.txt',header=T)
NewCode = UKB_crosswalk$NewCode[match(UKB_row_idx, as.numeric(UKB_crosswalk$NewCode_num)-100000000)]
UKB_embeddings = cbind(NewCode,UKB_embeddings)
UKB_embeddings$NewCode = sub('Org_','',UKB_embeddings$NewCode)

#### MGI embedding ####
MGI_row_idx_loc = paste0("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/",MGI_embedding_file,"/CodeEmbedding/New_",MGI_optimalWindow,"days_",MGI_optimalVecLen,"vec.RData")
load(MGI_row_idx_loc)
MGI_row_idx = as.numeric(row.names(CodeEmbedVec))
MGI_embeddings = read_csv(paste0("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/",MGI_embedding_file,"/CodeEmbedding/New_",MGI_optimalWindow,"days_",MGI_optimalVecLen,"vec.csv"))
MGI_crosswalk = read.table(paste0("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/",MGI_embedding_file,"/CodeCrosswalk/NewCode_num_10.txt"),header=T)
NewCode = MGI_crosswalk$NewCode[match(MGI_row_idx, as.numeric(MGI_crosswalk$NewCode_num)-100000000)]
MGI_embeddings = cbind(NewCode,MGI_embeddings)
MGI_embeddings$NewCode = sub('Org_','',MGI_embeddings$NewCode)

UKB_embeddings = UKB_embeddings[,1:3]
MGI_embeddings = MGI_embeddings[,1:3]

phecode_all_df = read.csv("/net/snowwhite/home/jiacong/EHRmapping/mapping2/data/MGI_UKB_freq.csv",header = TRUE)
phecode_all_MGI = phecode_all_df$phecode_MGI
phecode_all_UKB = phecode_all_df$phecode_UKB
phecode_all = union(unique(phecode_all_MGI),unique(phecode_all_UKB))

auc.before.all = NULL

for(i in 2:length(phecode_all)){

    phecode = phecode_all[i]

    UKB_phecode_ICD10 = find_ICD(phecode=phecode, mappingFile=icd_phecode,poolICD = UKB_embeddings$NewCode, ICDversion="ICD10")
    MGI_phecode_ICD10 = find_ICD(phecode=phecode, mappingFile=icdcm_phecode,poolICD = MGI_embeddings$NewCode, ICDversion = "ICD10")

    if(length(MGI_phecode_ICD10)==0|length(UKB_phecode_ICD10)==0){
        print("no ICD code found")
        next
    }

    # create the design matrix
    UKB_ids = unique(UKB_ICD$Deid_ID_num)
    MGI_ids = unique(MGI_ICD$Deid_ID_num)
    # UKB
    UKB_data = NULL
    for(j in 1:length(UKB_phecode_ICD10)){
    code = UKB_phecode_ICD10[j]
    temp = UKB_ICD[which(UKB_ICD$ICD==code),]
    UKB_temp = as.numeric(UKB_ids %in% temp$Deid_ID_num)
    UKB_data = cbind(UKB_data, UKB_temp)
    }
    colnames(UKB_data) = UKB_phecode_ICD10
    # MGI
    MGI_data = NULL
    for(j in 1:length(MGI_phecode_ICD10)){
    code = MGI_phecode_ICD10[j]
    temp = MGI_ICD[which(MGI_ICD$ICD==code),]
    MGI_temp = as.numeric(MGI_ids %in% temp$Deid_ID_num)
    MGI_data = cbind(MGI_data, MGI_temp)
    }
    colnames(MGI_data) = MGI_phecode_ICD10

    #### orginal cosine value ####
    codes_overlapped = intersect(UKB_phecode_ICD10,MGI_phecode_ICD10)
    if(length(codes_overlapped)==0){
        print("no overlapped ICD code found")
        next
    }
    MGI_data_overlapped = MGI_data[,codes_overlapped,drop=FALSE]
    UKB_data_overlapped = UKB_data[,codes_overlapped,drop=FALSE]
    data_combined_overlapped = as.data.frame(rbind(cbind(y=1,MGI_data_overlapped),cbind(y=0,UKB_data_overlapped)))
    # AUC
    model.formula = as.formula(paste0("y~",paste0(codes_overlapped,collapse = "+")))
    model = glm(model.formula,data = data_combined_overlapped,family="binomial")
    predicted <- predict(model, data_combined_overlapped, type="response")
    auc.before.est = round(as.numeric(ci.auc(data_combined_overlapped[,1], predicted)),3)
    (auc.before = paste0(auc.before.est[2]," (", auc.before.est[1],", ",auc.before.est[3], ") "))
    auc.before.i = c(phecode,auc.before.est[2],auc.before.est[1],auc.before.est[3])
    
    print(i)
    print(phecode)
    print(auc.before)
    auc.before.all = rbind(auc.before.all,auc.before.i)
}

write.csv(auc.before.all,paste0(outpath,"Step4_1_1_auc_before_all.csv"),row.names = FALSE)
