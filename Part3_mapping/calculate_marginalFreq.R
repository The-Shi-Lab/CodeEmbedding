# MGI, without repeated measurements
MGI_freqFile = fread("~/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis_noDuplicates/MGI_NewCode_freq_20221018_10.txt")
MGI_crosswalk = read.table('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis_noDuplicates/CodeCrosswalk/NewCode_num_10.txt',header=T)
MGI_freqFile$NewCode = MGI_crosswalk$NewCode[match(MGI_freqFile$NewCode_num, MGI_crosswalk$NewCode_num)]
MGI_freqFile$icd = sub(MGI_freqFile$NewCode, pattern = "Org_", replacement = "")
write.csv(MGI_freqFile, "~/EHRmapping/mapping2/data/MGI_mainDiagnosis_noDuplicates_freq.csv", row.names = FALSE)


# MGI, with repeated measurements
MGI_freqFile = fread("~/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/MGI_NewCode_freq_20221018_10.txt")
MGI_crosswalk = read.table('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/CodeCrosswalk/NewCode_num_10.txt',header=T)
MGI_freqFile$NewCode = MGI_crosswalk$NewCode[match(MGI_freqFile$NewCode_num, MGI_crosswalk$NewCode_num)]
MGI_freqFile$icd = sub(MGI_freqFile$NewCode, pattern = "Org_", replacement = "")
write.csv(MGI_freqFile, "~/EHRmapping/mapping2/data/MGI_mainDiagnosis_freq.csv", row.names = FALSE)

