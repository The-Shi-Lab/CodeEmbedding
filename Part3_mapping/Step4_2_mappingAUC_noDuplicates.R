
library(pROC)

load("~/EHRmapping/mapping2/data/designMatrix_org.RData")
dim(UKB_ICD)
dim(MGI_ICD)


#### functions ####
mappingAUC = function(cosine_matrix,source_designmatrix,target_designMatrix){
    Gamma_hat = apply(cosine_matrix,2,as.numeric)
    MGI_trans = as.matrix(MGI_data) %*% Gamma_hat
    data_comb = as.data.frame(rbind(cbind(y=1,MGI_trans),cbind(y=0,UKB_data)))
    colnames(data_comb) = c("y",UKB_phecode_ICD10)
    # AUC
    model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
    model = glm(model.formula,data = data_comb,family="binomial")
    summary(model)
    predicted <- predict(model, data_comb, type="response")
    auc.est = round(as.numeric(ci.auc(data_comb[,1], predicted)),3)
    result = auc = paste0(auc.est[2]," (", auc.est[1],", ",auc.est[3], ") ")
    return(list("auc.est" = auc.est,
                "result" = result))
}

get_cosine_matrix = function(sparsify,orgCos,truncatedCos){
    if (sparsify=="noSparsify"){
        cosine_matrix = orgCos
    } else if(sparsify == "top1"){
        cosine_matrix = t(apply(orgCos,1,function(x){x[x!=max(x)] = 0;x}))
    } else if (sparsify=="dataDrivenThreshold"){
        cosine_matrix = truncatedCos
    }
    return(cosine_matrix)
}

########################
# read the cosine matrix
MGI_embed_duplicates = FALSE
phecode = 401
outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"
cosineType = "orgCosine" 
cosineType = "refinedCosine"
# sparsify method 
sparsify = "top1"
# sparsify = "noSparsify"
# sparsify = "dataDrivenThreshold"

orgCos = read.csv(paste0("~/EHRmapping/mapping2/results/phecode_401_",cosineType,"_noDuplicates.csv"),header = TRUE, row.names = 1)
truncatedCos = read.csv(paste0("~/EHRmapping/mapping2/results/phecode_401_",cosineType,"_noDuplicates_thres.csv"),header = TRUE, row.names = 1)

UKB_phecode_ICD10 = colnames(orgCos)
MGI_phecode_ICD10 = rownames(orgCos)

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

(zeroEntry_MGI_all = mean(apply(MGI_data,1,sum)==0))
(zerEntry_UKB_all = mean(apply(UKB_data,1,sum)==0))


#### orginal cosine value ####
codes_overlapped = intersect(UKB_phecode_ICD10,MGI_phecode_ICD10)
MGI_data_overlapped = MGI_data[,codes_overlapped]
UKB_data_overlapped = UKB_data[,codes_overlapped]
data_combined_overlapped = as.data.frame(rbind(cbind(y=1,MGI_data_overlapped),cbind(y=0,UKB_data_overlapped)))
# AUC
model.formula = as.formula(paste0("y~",paste0(codes_overlapped,collapse = "+")))
model = glm(model.formula,data = data_combined_overlapped,family="binomial")
summary(model)
predicted <- predict(model, data_combined_overlapped, type="response")
head(predicted);head(data_combined_overlapped[,1]);head(data_combined_overlapped)
auc.before.est = round(as.numeric(ci.auc(data_combined_overlapped[,1], predicted)),3)
(auc.before = paste0(auc.before.est[2]," (", auc.before.est[1],", ",auc.before.est[3], ") "))

###### after harmonization ######
# cos_noSparsify = get_cosine_matrix(sparsify="noSparsify",orgCos,truncatedCos)
# (auc_noSparsify = mappingAUC(cos_noSparsify,source_designmatrix = MGI_data, target_designMatrix = UKB_data))
cos_top1 = get_cosine_matrix(sparsify="top1",orgCos,truncatedCos)
(auc_top1 = mappingAUC(cos_top1,source_designmatrix = MGI_data, target_designMatrix = UKB_data))
cos_dataDrivenThreshold = get_cosine_matrix(sparsify="dataDrivenThreshold",orgCos,truncatedCos)
(auc_dataDrivenThreshold = mappingAUC(cos_dataDrivenThreshold,source_designmatrix = MGI_data, target_designMatrix = UKB_data))

auc.before
auc_top1
auc_dataDrivenThreshold

result = data.frame("before_H_mean" = auc.before.est[2],"before_H_low" = auc.before.est[1],"before_H_up" = auc.before.est[3],
                    "top1_H_mean" = auc_top1$auc.est[2],"top1_H_low" =  auc_top1$auc.est[1],"top1_H_up" =  auc_top1$auc.est[3],
                    "dataDrivenThreshold_H_mean" = auc_dataDrivenThreshold$auc.est[2],"dataDrivenThreshold_H_low" =  auc_dataDrivenThreshold$auc.est[1],"dataDrivenThreshold_H_up" =  auc_dataDrivenThreshold$auc.est[3])
write.csv(result,paste0("~/EHRmapping/mapping2/results/phecode_",phecode,"_",cosineType,"_noDuplicates_AUC.csv"),row.names = FALSE)

