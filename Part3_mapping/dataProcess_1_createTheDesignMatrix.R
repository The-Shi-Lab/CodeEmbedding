library(data.table)
library(dplyr)

outfile_Rdata = "~/EHRmapping/mapping2/data/designMatrix_org.RData"

########### MGI_data ##########
MGI_NewCode_time_org = fread("~/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis_noDuplicates/MGI_NewCode_time_20221018_10.txt")
MGI_NewCode_crosswalk = fread("~/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis_noDuplicates/CodeCrosswalk/NewCode_num_source_10.txt")

MGI_NewCode_time_1 = MGI_NewCode_time_org[!duplicated(MGI_NewCode_time_org),]
MGI_NewCode_time_2 = MGI_NewCode_time_org[!duplicated(MGI_NewCode_time_org[,c("Deid_ID_num","NewCode_num")]),]
dim(MGI_NewCode_time_org)
dim(MGI_NewCode_time_1)
dim(MGI_NewCode_time_2)

MGI_ICD = MGI_NewCode_time_2 %>%
                     mutate(prefixICD = MGI_NewCode_crosswalk$NewCode[match(NewCode_num, MGI_NewCode_crosswalk$NewCode_num)])%>%
                     mutate(ICD = sub('Org_','',prefixICD))%>%
                     select(Deid_ID_num,ICD)
############ UKB data ##########
UKB_NewCode_time_org = fread("~/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/UKB_NewCode_time_20220927_5.txt")
UKB_NewCode_crosswalk = fread("~/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeCrosswalk/NewCode_num_source_5.txt")

UKB_NewCode_time_1 = UKB_NewCode_time_org[!duplicated(UKB_NewCode_time_org),]
UKB_NewCode_time_2 = UKB_NewCode_time_org[!duplicated(UKB_NewCode_time_org[,c("Deid_ID_num","NewCode_num")]),]
dim(UKB_NewCode_time_org)
dim(UKB_NewCode_time_1)
dim(UKB_NewCode_time_2)

UKB_ICD = UKB_NewCode_time_2 %>%
                     mutate(prefixICD = UKB_NewCode_crosswalk$NewCode[match(NewCode_num, UKB_NewCode_crosswalk$NewCode_num)])%>%
                     mutate(ICD = sub('Org_','',prefixICD))%>%
                     select(Deid_ID_num,ICD)


save(MGI_ICD, UKB_ICD,file = outfile_Rdata)
