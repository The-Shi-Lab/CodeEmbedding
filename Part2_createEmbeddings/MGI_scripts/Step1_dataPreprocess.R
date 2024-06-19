library(data.table)
library(dplyr)
library(stringr)


source("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/scripts/function.harmonizeICD9.r")
threshold = 10
path = "/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/"
today = "20221018"

######
# dx_org <- fread('/net/junglebook/magic_data/ehr_data_20210318/HPI5635_Diagnosis_DEIDENTIFIED.txt', select=c("DeId_PatientID","DaysSinceBirth","DiagnosisCode","Lexicon"))

dx_org0 <- fread('/net/junglebook/magic_data/ehr_data_20210318/HPI5635_Diagnosis_DEIDENTIFIED.txt')
dx_org0 = dx_org0[,c("DeId_PatientID","DaysSinceBirth","DiagnosisCode","Lexicon","DiagnosisCodeType","DiagnosisSource")]
dim(dx_org0) # 117296568

# dx_org0[grepl("I63",dx_org0$DiagnosisCode),]


dx_org0 = dx_org0 %>%
    filter(DiagnosisCodeType=="Primary"|DiagnosisCodeType=="Secondary"|DiagnosisSource=="MiChart Visit Diagnosis") %>%
    select(-DiagnosisCodeType,-DiagnosisSource)
dim(dx_org0) # 53549849, 45%

# note on armis2, the column name for patient id is Deid_ID
dx = dx_org0
colnames(dx) <- c('Deid_ID', 'DaysSinceBirth', 'DxCode',"Lexicon")

# remove missing date
dx$DaysSinceBirth = as.integer(dx$DaysSinceBirth)
dx = na.omit(dx)
dim(dx) # 53396917        4

# remove duplicates
dx = dx[!duplicated(dx), ]
dim(dx) # 34413035        4

# remove IMO0001 and IMO0002
dx = dx[(dx$DxCode!= "IMO0001") & (dx$DxCode!= "IMO0002"),]
dim(dx) # 34393713        4

dx = dx[,c("Deid_ID","DaysSinceBirth","DxCode")]
dx = dx %>%
  dplyr::group_by(Deid_ID) %>% 
  dplyr::arrange(DaysSinceBirth,.by_group = TRUE)
# ##### remove repeated measurements #####
# dx_person = dx[,c("Deid_ID","DxCode")]
# duplicates_ID = which(duplicated(dx_person))
# dx = dx[-duplicates_ID,]
# print("The dimension of no-duplicated data is")
# print(dim(dx))

####################################
######## handle rare codes #########
####################################

####  source phecode files ####
#### source matched icd-phecodes files ######
icd9_phecode <- read.csv('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/data/phecode_icd9_map_unrolled.csv')
icd10_phecode <- read.csv('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/data/Phecode_map_v1_2_icd10cm_beta.csv')
colnames(icd9_phecode) <- c('ICD_dx', 'phecode')
icd10_phecode <- icd10_phecode[,c('icd10cm','phecode')]
colnames(icd10_phecode) <- c('ICD_dx', 'phecode')
icd_phecode <- rbind(icd9_phecode, icd10_phecode)


#### source the UKB tree structure files ######
ICD9codes <- fread("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/data/coding87.tsv")
ICD9codes[!grepl("Block",coding),ICD:=sapply(coding,harmonizeICD9)]
ICD10codes <- fread("/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/data/coding19.tsv")
ICD10codes <- ICD10codes[!grepl("Block",coding),]
ICD10codes[,ICD10category:=gsub("([A-Z][0-9]{2}).+","\\1",coding)]
ICD10codes[,ICD10suffix:=gsub("[A-Z].+","",gsub("^[A-Z][0-9]{2}","",coding))]
ICD10codes[,ICD:=paste0(ICD10category,ifelse(ICD10suffix == "","","."),ICD10suffix)]
ICD10codes <- ICD10codes[,c("ICD10category","ICD10suffix"):=NULL]
ICD_Dx = rbind(ICD9codes,ICD10codes)

#### start #####
dx = dx[,c("Deid_ID","DaysSinceBirth","DxCode")]
dx_freq = as.data.frame(table(dx$DxCode))
colnames(dx_freq)[1] = "DxCode"
CodeAllInfo=merge(dx,dx_freq,by="DxCode")
CodeAllInfo$NewCode=''
CodeAllInfo$Source=''

FrequentCode=CodeAllInfo[CodeAllInfo$Freq>=threshold,]
FrequentCode$NewCode=paste0("Org_",FrequentCode$DxCode)
FrequentCode$Source="Original"

RareCode=CodeAllInfo[!CodeAllInfo$Freq>=threshold,]

# digit = 1
RareCode1 = RareCode %>% 
  dplyr::mutate(NewCode=icd_phecode$phecode[match(DxCode,icd_phecode$ICD_dx)]) %>%
  dplyr::mutate(Source=ifelse(!is.na(NewCode), 'Phecode_org', '')) %>%
  dplyr::mutate(NewCode1= if_else((Source=="")&!is.na(str_extract(DxCode,".*\\.[0-9]")),str_extract(DxCode,".*\\.[0-9]"), as.character(NewCode))) %>%
  dplyr::mutate(Source1 = if_else((Source=="")&!is.na(str_extract(DxCode,".*\\.[0-9]")), "ICD_trc", Source)) %>% 
  dplyr::mutate(NewCode2 = if_else(Source1=="", as.character(DxCode), NewCode1)) %>%
  dplyr::mutate(Source2 = if_else(Source1=="", "ICD_org", Source1)) %>%
  dplyr::mutate(NewCode3 = ifelse(Source2=="Phecode_org", paste0("Porg_",NewCode2), ifelse(
    Source2=="ICD_trc", paste0("Itrc_",NewCode2), paste("Iorg_",NewCode2)  
  ))) %>%
  dplyr::select("Deid_ID","DaysSinceBirth","DxCode","Freq","NewCode2","NewCode3","Source2")
colnames(RareCode1)[5:7] = c("NewCode","NewCode_prefix","Source")
#### note: since the second re-coding needs information in NewCode, we keep it. ####

# re-calculate the frequency
RareCode2 = RareCode1
RareCode2_freq = as.data.frame(table(RareCode2$NewCode_prefix))
colnames(RareCode2_freq)[1] = "NewCode"
RareCode2$Freq2 = RareCode2_freq$Freq[match(RareCode2$NewCode_prefix,RareCode2_freq$NewCode)]
table(RareCode2$Source[which(RareCode2$Freq2<threshold)])

RareCode3 = RareCode2 %>%
  dplyr::mutate(NewCode1 = if_else((Freq2<threshold)&(Source=="Phecode_org")&(!is.na(str_extract(NewCode,".*(?=\\.)"))),
                                   str_extract(NewCode,".*(?=\\.)"), NewCode) ) %>%
  dplyr::mutate(Source1 = if_else((Freq2<threshold)&(Source=="Phecode_org")&(!is.na(str_extract(NewCode,".*(?=\\.)"))),
                                  "Phecode_trc", Source) ) %>%
  dplyr::mutate(NewCode_rollup = if_else( (Freq2<threshold)&(Source=="ICD_trc"), str_extract(NewCode,".*(?=\\.)"), NewCode1))%>%
  dplyr::mutate(NewCode2 = if_else((Freq2<threshold)&(Source%in%c("ICD_trc","ICD_org")), 
                                   as.character(ICD_Dx$parent_id[match(NewCode1,ICD_Dx$ICD)]), as.character(NewCode1))) %>%
  dplyr::mutate(Source2 =  if_else((Freq2<threshold)&(Source%in%c("ICD_trc","ICD_org")), 
                                   "ICD_group",Source1)) %>%
  dplyr::mutate(NewCode3 = ifelse(Source2=="Phecode_org", paste0("Porg_",NewCode2), ifelse(
    Source2=="Phecode_trc", paste0("Ptrc_",NewCode2), ifelse(
      Source2=="ICD_trc", paste0("Itrc_",NewCode2), paste0("Igrp_",NewCode2)
    ) ))) %>%
  dplyr::select("Deid_ID","DaysSinceBirth","DxCode","NewCode3","Source2")
colnames(RareCode3)[4:5] = c("NewCode","Source")

RareCode3_freq = as.data.frame(table(RareCode3$NewCode))
colnames(RareCode3_freq)[1] = "NewCode"
RareCode3$Freq = RareCode3_freq$Freq[match(RareCode3$NewCode,RareCode3_freq$NewCode)]
table(RareCode3$Source[which(RareCode3$Freq<threshold)])

# combine frequent codes and rare codes
NewCodeInfo=rbind(FrequentCode,RareCode3)
NewCodeInfo=NewCodeInfo[,c("Deid_ID","DaysSinceBirth", "Freq","NewCode","Source")]
NewCode_final = NewCodeInfo[-which(NewCodeInfo$Freq<threshold),]

##### check deleted codes ######
# new_code_below=NewCodeInfo[which(NewCodeInfo$Freq<threshold),]
# new_code_below=new_code_below[!duplicated(new_code_below$NewCode),]
# (proportion=nrow(new_code_below)/length(unique(NewCodeInfo$NewCode)) ) # proportion of codes<threshold
# 
# frequency_table = new_code_below %>%
#   group_by(Freq,Source) %>%
#   summarize(Source_freq=n())
# 
# ggplot(frequency_table, aes(fill=Source, y=Source_freq, x=factor(Freq),label = Source_freq)) + 
#   geom_bar(position="stack", stat="identity")+
#   geom_text(size = 5, position = position_stack(vjust = 0.5))+
#   xlab("Frequency")+
#   ylab("Number of NewCode")+
#   ggtitle(paste0("UKB: Threshold=",threshold,", Truncated digits=", digit))+
#   theme_bw()+  
#   annotate(geom="text", label=paste0(round(proportion*100,1),"%"),
#            x = -Inf, y = Inf, hjust = -0.5, vjust = 1.2,size=10, color="red")

####### Re-number patient id and NewCodes ######

# rm(list=setdiff(ls(),"NewCode_final"))
# threshold = 10
# path = "/nfs/turbo/mgi-shixu/project/MGI_embedding/MGI_thres10/"
# today = "20221018"

new_unique = unique(NewCode_final$NewCode)
NewCode_num <- data.frame(NewCode=as.character(new_unique), 
                          NewCode_num=as.character((1:length(new_unique)+100000000)))
NewCode_final$NewCode_num = NewCode_num$NewCode_num[match(NewCode_final$NewCode, NewCode_num$NewCode)]

Deid_ID_unique = unique(NewCode_final$Deid_ID)
Deid_ID_num = data.frame(Deid_ID=as.character(Deid_ID_unique), 
                         Deid_ID_num=as.character((1:length(Deid_ID_unique))))
NewCode_final$Deid_ID_num = Deid_ID_num$Deid_ID_num[match(NewCode_final$Deid_ID, Deid_ID_num$Deid_ID)]

####### outputs ##############
# output 1: 2 columns: NewCode, NewCode_num
write.table(NewCode_num, file=paste0(path,'CodeCrosswalk/NewCode_num_',threshold,'.txt'), row.names=F)

# output 2: 3 columns: NewCode, NewCode_num, Source
NewCode_num_Source = NewCode_num %>%
  mutate(Source = NewCode_final$Source[match(NewCode,NewCode_final$NewCode)])
table(NewCode_num_Source$Source)
write.table(NewCode_num_Source, file=paste0(path,'CodeCrosswalk/NewCode_num_source_',threshold,'.txt'), row.names=F)
# NewCode_num_Source = fread(paste0(path,'CodeCrosswalk/NewCode_num_source_',threshold,'.txt'))
# table(NewCode_num_Source$Source)
# NewCode_num_Source = fread(paste0('~/EHRmapping/createUKBphenome/UKB_icdcode_thres5_digit1_pmi/CodeCrosswalk/NewCode_num_source_',5,'.txt'))
# table(NewCode_num_Source$Source)

# output 3: 2 columns: NewCode_num, Freq
new_freq=as.data.frame(table(NewCode_final$NewCode_num))
colnames(new_freq)=c("NewCode_num","Freq")
sum(new_freq$Freq<threshold) # check the frequency
fwrite(new_freq,paste0(path,"MGI_NewCode_freq_",today,'_',threshold,".txt"),sep="\t",row.names=F,col.names=T,quote=T)

rm(list=setdiff(ls(),"NewCode_final"))
threshold = 10
path = "/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/"
today = "20221018"

# output 4: 3 columns: Deid_ID_num,DaysSinceBirth, NewCode_num
NewCode_output<-NewCode_final %>% 
  dplyr::select(c("Deid_ID_num","DaysSinceBirth", "NewCode_num")) %>%
  dplyr::group_by(Deid_ID_num) %>% 
  dplyr::arrange(DaysSinceBirth,.by_group = TRUE)
colnames(NewCode_output) = c("Deid_ID_num","DaysSinceBirth", "NewCode_num")
fwrite(NewCode_output,paste0(path,"MGI_NewCode_time_",today,'_',threshold,".txt"),sep="\t",row.names=F,col.names=T,quote=T)

