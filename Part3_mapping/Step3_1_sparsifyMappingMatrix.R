
MGI_embed_duplicates = TRUE
phecode = 401
outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"
cosineType = "orgCosine" # refinedCosine
cosineType = "refinedCosine"

source("~/EHRmapping/mapping2/sourceFun/function_mapping.R")
source("~/EHRmapping/mapping2/sourceFun/refine_mapping.R")

load(paste0(outpath,"phecode_",phecode,"_orgCosine_withDuplicates.RData"))
# load(paste0(outpath,"phecode_",phecode,"_refinedCosine_withDuplicates_freq.RData"))


coss = read.csv(paste0("/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/phecode_",phecode,"_",cosineType,"_withDuplicates.csv"),header = TRUE, row.names = 1)

# top-1 thresholding
coss_thres_top1 = t(apply(coss,1,function(x){ (x==max(x))*x }))

if(MGI_embed_duplicates){
    outfile = paste0(outpath,"phecode_",phecode,"_",cosineType,"_withDuplicates_top1.csv")
} else {
    outfile = paste0(outpath,"phecode_",phecode,"_",cosineType,"_noDuplicates_top1.csv")}
write.csv(coss_thres_top1, outfile, row.names = TRUE)

# data-drive thresholding
thresholds = seq(0.1,0.7,by=0.01)
MGI_codes = rownames(coss)
UKB_codes = colnames(coss)
MGI_embeddings_aligned = MGI_embeddings_map %*% beta_hat

if (cosineType == "orgCosine"){
    loss.evaluate = coss_cv(y_hat = MGI_embeddings_aligned,
                            y = UKB_embeddings_map,
                            thres = thresholds,
                            marginal=FALSE, UKB_freq_norm=MGI_freq, MGI_freq_norm=UKB_freq,
                            nfold=20,
                            leaveOneOut=FALSE)
}else if(cosineType=="refinedCosine"){
    loss.evaluate = coss_cv(y_hat = MGI_embeddings_aligned,
                            y = UKB_embeddings_map,
                            thres = thresholds,
                            marginal=TRUE, UKB_freq_norm=MGI_freq, MGI_freq_norm=UKB_freq,
                            nfold=20,
                            leaveOneOut=FALSE)
}

(loss.results = loss.evaluate$loss)
(cutoff.ICD10_SR = loss.evaluate$cutoff.min)

coss_truncated = coss
coss_truncated[which(coss_truncated<cutoff.ICD10_SR,arr.ind = TRUE)]=0

if(MGI_embed_duplicates){
    outfile = paste0(outpath,"phecode_",phecode,"_",cosineType,"_withDuplicates_thres.csv")
} else {
    outfile = paste0(outpath,"phecode_",phecode,"_",cosineType,"_noDuplicates_thres.csv")}
write.csv(coss_truncated, outfile, row.names = TRUE)
