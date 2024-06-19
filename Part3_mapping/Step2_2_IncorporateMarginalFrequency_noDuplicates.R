
source("~/EHRmapping/mapping2/sourceFun/function_mapping.R")
source("~/EHRmapping/mapping2/sourceFun/refine_mapping.R")

outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"
MGI_embed_duplicates = FALSE
MGI_freqFile = read.csv("~/EHRmapping/mapping2/data/MGI_mainDiagnosis_noDuplicates_freq.csv")[,-1]
UKB_freqFile = read.csv("~/EHRmapping/mapping2/data/UKB_ICD10_noDuplicates_freq.csv")[,-1]

# read coss file
coss = read.csv("/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/phecode_401_orgCosine_noDuplicates.csv",header = TRUE, row.names = 1)

# get the frequency
MGI_codes = rownames(coss)
UKB_codes = colnames(coss)
MGI_freq = MGI_freqFile[match(MGI_codes,MGI_freqFile$icd),"freq"]
MGI_freq = MGI_freq/sum(MGI_freq)
UKB_freq = UKB_freqFile[match(UKB_codes,UKB_freqFile$icd),"freq"]
UKB_freq = UKB_freq/sum(UKB_freq)

# update the cosine matrix
coss_refined = refine_map_mt(map_mt = as.matrix(coss), source = MGI_freq, target = UKB_freq)

if(MGI_embed_duplicates){
    outfile = paste0(outpath,"phecode_",phecode,"_refinedCosine_withDuplicates.csv")
} else {
    outfile = paste0(outpath,"phecode_",phecode,"_refinedCosine_noDuplicates.csv")}

write.csv(coss_refined, outfile, row.names = TRUE)

