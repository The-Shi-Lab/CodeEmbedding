MGI_freqFile = read.csv("~/EHRmapping/mapping2/data/MGI_mainDiagnosis_noDuplicates_freq.csv")[,-1]
UKB_freqFile = read.csv("~/EHRmapping/mapping2/data/UKB_ICD10_noDuplicates_freq.csv")[,-1]

icd_phecode = read.csv("~/EHRmapping/mapping/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping/data/icdcm_phecodes.csv")[,-1]

MGI_freqFile0 = MGI_freqFile %>%
               mutate(S=ifelse(grepl("Org_",NewCode),1,0)) %>%
               mutate(icd = sub("Org_","",NewCode)) %>%
               filter(S==1) %>%
               mutate(phecode = floor(icdcm_phecode[match(icd,icdcm_phecode$ICD_dx),"phecode"])) %>%
               select(icd,phecode,freq)
head(MGI_freqFile0)

UKB_freqFile0 = UKB_freqFile %>%
                mutate(phecode = floor(icd_phecode[match(icd,icd_phecode$ICD_dx),"phecode"])) %>%
                select(icd,phecode,freq)

data_merged = merge(MGI_freqFile0,UKB_freqFile0,by="icd",all=TRUE)
data_merged1 = data_merged
colnames(data_merged1) = c("icd","phecode_MGI","freq_MGI","phecode_UKB","freq_UKB")
head(data_merged1)


n_MGI = 80000
n_UKB = 500000
data_merged2 = data_merged1 %>%
               filter(grepl("^[A-Z]",icd)) %>%
               group_by(phecode_MGI) %>%
               mutate(MGI_prevalence = sum(freq_MGI,na.rm = TRUE)/n_MGI,
                      UKB_prevalence = sum(freq_UKB,na.rm = TRUE)/n_UKB,  
                      freq_MGI_norm = freq_MGI/n_MGI,
                      freq_UKB_norm = freq_UKB/n_UKB) %>%
            #    filter(!is.na(phecode_MGI)&!is.na(phecode_UKB)) %>%
               filter(!is.na(phecode_UKB)) %>%
               arrange(phecode_UKB) 
data_merged2
data_merged2[,c("icd","phecode_UKB","MGI_prevalence","UKB_prevalence","freq_MGI_norm","freq_UKB_norm")]

write.csv(data_merged2,"~/EHRmapping/mapping2/data/MGI_UKB_freq.csv",row.names = FALSE)
