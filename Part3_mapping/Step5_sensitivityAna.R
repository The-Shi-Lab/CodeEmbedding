source("~/EHRmapping/mapping2/sourceFun/function_mapping.R")
source("~/EHRmapping/mapping2/sourceFun/refine_mapping.R")


# Note: sensitivity analysis is to incorporate marginal frequency with no duplicates
load(paste0(outpath,"phecode_",phecode,"_orgCosine_withDuplicates.RData"))

############ repeat step2: incorporate marginal frequency ############
outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"
MGI_embed_duplicates = TRUE
phecode = 401
MGI_freqFile = read.csv("~/EHRmapping/mapping2/data/MGI_mainDiagnosis_noDuplicates_freq.csv")[,-1]
# MGI_freqFile = read.csv("~/EHRmapping/mapping2/data/MGI_mainDiagnosis_freq.csv")[,-1]
head(MGI_freqFile); colnames(MGI_freqFile) = c("freq","NewCode","icd")
UKB_freqFile = read.csv("~/EHRmapping/mapping2/data/UKB_ICD10_noDuplicates_freq.csv")[,-1]

# read coss file
coss = read.csv("/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/phecode_401_orgCosine_withDuplicates.csv",header = TRUE, row.names = 1)

# get the frequency
MGI_codes = rownames(coss)
UKB_codes = colnames(coss)
MGI_freq = MGI_freqFile[match(MGI_codes,MGI_freqFile$icd),"freq"]
MGI_freq = MGI_freq/sum(MGI_freq)
UKB_freq = UKB_freqFile[match(UKB_codes,UKB_freqFile$icd),"freq"]
UKB_freq = UKB_freq/sum(UKB_freq)

# update the cosine matrix
coss_refined = refine_map_mt(map_mt = as.matrix(coss), source = MGI_freq, target = UKB_freq)

############ repeat step3: sparsitfy the cosine matrix ############
cosineType = "refinedCosine"

source("~/EHRmapping/mapping2/sourceFun/function_mapping.R")
source("~/EHRmapping/mapping2/sourceFun/refine_mapping.R")

coss = coss_refined
# top-1 thresholding
coss_thres_top1 = t(apply(coss,1,function(x){ (x==max(x))*x }))
outfile = paste0(outpath,"sensitivityAna_phecode_",phecode,"_",cosineType,"_withDuplicates_top1.csv")
write.csv(coss_thres_top1, outfile, row.names = TRUE)

# data-drive thresholding
thresholds = seq(0.1,0.7,by=0.01)
MGI_codes = rownames(coss)
UKB_codes = colnames(coss)
MGI_embeddings_aligned = MGI_embeddings_map %*% beta_hat

loss.evaluate = coss_cv(y_hat = MGI_embeddings_aligned,
                        y = UKB_embeddings_map,
                        thres = thresholds,
                        marginal=TRUE, UKB_freq_norm=MGI_freq, MGI_freq_norm=UKB_freq,
                        nfold=20,
                        leaveOneOut=FALSE)
(loss.results = loss.evaluate$loss)
(cutoff.ICD10_SR = loss.evaluate$cutoff.min)

coss_truncated = coss
coss_truncated[which(coss_truncated<cutoff.ICD10_SR,arr.ind = TRUE)]=0
outfile = paste0(outpath,"sensitivityAna_phecode_",phecode,"_",cosineType,"_withDuplicates_thres.csv")
write.csv(coss_truncated, outfile, row.names = TRUE)

######################### repeat step4: calculate the AUC #########################

library(pROC)
load("~/EHRmapping/mapping2/data/designMatrix_org.RData")

#### functions ####
mappingAUC = function(cosine_matrix,source_designmatrix,target_designMatrix){
    Gamma_hat = apply(cosine_matrix,2,as.numeric)
    MGI_trans = as.matrix(MGI_data) %*% Gamma_hat
    data_comb = as.data.frame(rbind(cbind(y=1,MGI_trans),cbind(y=0,UKB_data)))
    colnames(data_comb) = c("y",UKB_phecode_ICD10)
    # AUC
    model.formula = as.formula(paste0("y~",paste0(UKB_phecode_ICD10,collapse = "+")))
    model = glm(model.formula,data = data_comb,family="binomial")
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

orgCos = coss_refined
truncatedCos = coss_truncated

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


# #### orginal cosine value ####
# codes_overlapped = intersect(UKB_phecode_ICD10,MGI_phecode_ICD10)
# MGI_data_overlapped = MGI_data[,codes_overlapped]
# UKB_data_overlapped = UKB_data[,codes_overlapped]
# data_combined_overlapped = as.data.frame(rbind(cbind(y=1,MGI_data_overlapped),cbind(y=0,UKB_data_overlapped)))
# # AUC
# model.formula = as.formula(paste0("y~",paste0(codes_overlapped,collapse = "+")))
# model = glm(model.formula,data = data_combined_overlapped,family="binomial")
# predicted <- predict(model, data_combined_overlapped, type="response")
# auc.before.est = round(as.numeric(ci.auc(data_combined_overlapped[,1], predicted)),3)
# (auc.before = paste0(auc.before.est[2]," (", auc.before.est[1],", ",auc.before.est[3], ") "))



###### after harmonization ######
# cos_noSparsify = get_cosine_matrix(sparsify="noSparsify",orgCos,truncatedCos)
# (auc_noSparsify = mappingAUC(cos_noSparsify,source_designmatrix = MGI_data, target_designMatrix = UKB_data))
cos_top1 = get_cosine_matrix(sparsify="top1",orgCos,truncatedCos)
(auc_top1 = mappingAUC(cos_top1,source_designmatrix = MGI_data, target_designMatrix = UKB_data))
cos_dataDrivenThreshold = get_cosine_matrix(sparsify="dataDrivenThreshold",orgCos,truncatedCos)
(auc_dataDrivenThreshold = mappingAUC(cos_dataDrivenThreshold,source_designmatrix = MGI_data, target_designMatrix = UKB_data))

result = data.frame("top1_H_mean" = auc_top1$auc.est[2],"top1_H_low" =  auc_top1$auc.est[1],"top1_H_up" =  auc_top1$auc.est[3],
                    "dataDrivenThreshold_H_mean" = auc_dataDrivenThreshold$auc.est[2],"dataDrivenThreshold_H_low" =  auc_dataDrivenThreshold$auc.est[1],"dataDrivenThreshold_H_up" =  auc_dataDrivenThreshold$auc.est[3])

result

write.csv(result,paste0("~/EHRmapping/mapping2/results/sensitivityAna_phecode_",phecode,"_",cosineType,"_withDuplicates_AUC.csv"),row.names = FALSE)


############ draw the plot ############

cosineType = "refinedCosine"

orgCos = coss_refined
truncatedCos = coss_thruncated
title = "Sensitivity analysis: ICD10 mapping from MGI to UKB using the refined cosine"
file_name ="~/EHRmapping/mapping2/results/sensitivityAna_phecode_401_withDuplicates_refined.png"


UKB_phecode = read.csv('~/EHRmapping/mapping2/data/phecode_icd10.csv')
UKB_phecode = na.omit(UKB_phecode)
colnames(UKB_phecode)[1:3] = c("icd10","icd10.string","phecode")
MGI_phecode = read.csv('~/EHRmapping/mapping2/data/Phecode_map_v1_2_icd10cm_beta.csv')
colnames(MGI_phecode)[1:3] = c("icd10","icd10.string","phecode")

#### original cos ####
data = orgCos
cos = orgCos
MGI_codes = rownames(cos)
UKB_codes = colnames(cos)
cos = cbind(MGI_codes,cos)

data_truncated = truncatedCos
cos_truncated = as.data.frame(truncatedCos)
cos_truncated = cbind(MGI_codes, cos_truncated)


# set up the y-axis for all codes
all_codes = unique(c(UKB_codes,MGI_codes))
all_codes = all_codes[order(all_codes)]
codes_index = data.frame("index" = 1:length(all_codes),"codes"=all_codes)
codes_index$id = 1:nrow(codes_index)

cos_long = cos %>% pivot_longer(cols=colnames(cos)[-1],
                                names_to='UKB_codes',
                                values_to='cos')
cos_long$cos = round(as.numeric(cos_long$cos),3)

cos_truncated_long = cos_truncated %>% pivot_longer(cols=colnames(cos_truncated)[-1],
                                                    names_to='UKB_codes',
                                                    values_to='cos_tr')
cos_merge = merge(cos_long, cos_truncated_long, by=c("MGI_codes","UKB_codes"))
cos_merge2 = cos_merge %>%
  group_by(MGI_codes) %>%
  arrange(desc(cos),.by_group = TRUE)
cos_merge2$cos_order = rep(1:length(UKB_codes),length(MGI_codes))
cos_merge3 = cos_merge2 %>%
  mutate(top1 = ifelse(cos_order==1, cos, 0)) %>%
  mutate(top2 = ifelse(cos_order==2, cos, 0))
cos_merge3$cos_tr = as.numeric(cos_merge3$cos_tr)

cos_merge3$UKB_index = codes_index$id[match(cos_merge3$UKB_codes, codes_index$codes)]
cos_merge3$MGI_index = codes_index$id[match(cos_merge3$MGI_codes, codes_index$codes)]
cos_merge3$UKB_codeString = UKB_phecode$icd10.string[match(cos_merge3$UKB_codes, UKB_phecode$icd10)]
cos_merge3$MGI_codeString = MGI_phecode$icd10.string[match(cos_merge3$MGI_codes, MGI_phecode$icd10)]
cos_merge3$top2=0

color_data = data.frame("x" = c(1,1),'y'=c(1,1),"col"=c("Pass threshold","Top 1"))
(pic = ggplot(data=cos_merge3, aes(x=-1,y=MGI_index))+
  scale_x_continuous(limits=c(-24,14.5))+
  geom_point(data = color_data, aes(x=x,y=y,col=col))+
  scale_color_manual(values = c("blue","grey"),name="")+
  ggtitle(title)+
  guides(color = guide_legend(override.aes = list(size = 20)))+
  geom_segment(aes(x=-1, y=MGI_index, xend = 1, yend = UKB_index),color="grey",size = abs(cos_merge3$top1)*8)+ # top1
  geom_segment(aes(x=-1, y=MGI_index, xend = 1, yend = UKB_index),color="blue",linetype=3,size = abs(cos_merge3$cos_tr)*8)+ # truncated
  geom_point(size=8)+
  geom_text(data=cos_merge3, aes(y=UKB_index,label=UKB_codes),nudge_x = 3.5, size=12,hjust=1)+
  geom_text(data=cos_merge3, aes(y=UKB_index,label=UKB_codeString),nudge_x = 4, size=12,hjust = 0)+
  geom_point(data=cos_merge3,aes(x=1,y=UKB_index),size=8)+
  geom_text(data=cos_merge3, aes(y=MGI_index,label=MGI_codes),nudge_x = -1.1, size=12,hjust=0)+
  geom_text(data=cos_merge3, aes(y=MGI_index,label=MGI_codeString),nudge_x = -1.5, size=12,hjust = 1)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="bottom",
        legend.key=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.text = element_text(size=30),
        legend.key.size = unit(3, 'cm'),
        plot.title = element_text(size = 40, face = "bold"))
)



png(file_name, width = 4400, height = 1600)

print(pic)

dev.off()

