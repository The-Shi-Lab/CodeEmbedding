if(!file.exists("./UKB_icdcodes")) dir.create("./UKB_icdcode")
today <- format(Sys.Date(), format =  "%Y%m%d")

### functions from Lars to help preprocess UKB data ###
# extract information in different instances/arrays of each field and transform them into longformat 
source("./scripts/function.reformatUKB.r")
# format icd9 codes, like adding decimal points
source("./scripts/function.harmonizeICD9.r")
# format icd10 codes, like adding decimal points
source("./scripts/function.harmonizeICD10.r")

### obtain icd10 codes with time ###
ukb_icd10_diagnose_time=reformatUKB(fields=c(41270, 41280),dataCoding = T)
ukb_icd10_cancer_time=reformatUKB(fields=c(40006,40005),dataCoding = T)
ukb_icd10_death_primary_time=reformatUKB(fields=c(40001,40000),dataCoding = T)
# external icd10 codes and secondary cause of death icd10 codes do not have corresponding time
# but they will be useful if we consider patients' lifetime code pairs
ukb_icd10_external=reformatUKB(fields=41201,dataCoding = T)
ukb_icd10_death_second=reformatUKB(fields=40002,dataCoding = T)

# icd10 diagnostic codes
icd10_diagnose_time<-ukb_icd10_diagnose_time[[1]]
colnames(icd10_diagnose_time)<-c("id","icd10","time")
# delete quotation marks in the raw data
icd10_diagnose_time<-icd10_diagnose_time[,icd10:=gsub('\\"\\"(+.+)\\"\\"',"\\1",icd10)]
icd10_diagnose_time<-icd10_diagnose_time[,icd10:=sapply(icd10,harmonizeICD10)]

# icd10 cancer codes
icd10_cancer_time<-ukb_icd10_cancer_time[[1]]
colnames(icd10_cancer_time)<-c("id","icd10","time")
icd10_cancer_time[,icd10:=sapply(icd10,harmonizeICD10)]

# icd10 primary cause of death codes
icd10_death_primary_time<-ukb_icd10_death_primary_time[[1]]
colnames(icd10_death_primary_time)<-c("id","icd10","time")
icd10_death_primary_time[,icd10:=sapply(icd10,harmonizeICD10)]

# icd10 external diagnostic codes
icd10_external<-ukb_icd10_external[[1]]
colnames(icd10_external)<-c("id","icd10")
icd10_external[,icd10:=sapply(icd10,harmonizeICD10)]

# icd10 secondary cause of death codes
icd10_death_second<-ukb_icd10_death_second[[1]]
colnames(icd10_death_second)<-c("id","icd10")
icd10_death_second[,icd10:=sapply(icd10,harmonizeICD10)]

# combine icd10 codes from different sources
all_icd10_time<-dplyr::full_join(icd10_death_primary_time,icd10_cancer_time)
all_icd10_time<-dplyr::full_join(all_icd10_time,icd10_diagnose_time)

# make a NA-omitted copy for merging later
all_icd10_time_copy=na.omit(all_icd10_time)

# obtain time-free patientid-icd10code pairs data
all_icd10<-dplyr::full_join(all_icd10_time,icd10_external)
all_icd10<-dplyr::full_join(all_icd10,icd10_death_second)
all_icd10<-all_icd10[,c(1,2)]
dim(all_icd10)
# [1] 4240684       2

# remove missing value
all_icd10 <- all_icd10[!(is.na(all_icd10$id)|is.na(all_icd10$icd10)),]
all_icd10 <- all_icd10[order(all_icd10$id, all_icd10$icd10),]

# ensure each patientid-icd10code pair occur but only occur once
all_icd10<-all_icd10[!duplicated(all_icd10),]
dim(all_icd10)
# [1] 4187227       2

# add prefix "10_"for these icd10 codes to distinguish from icd9 codes
class(all_icd10$icd10)
#all_icd10$icd10<-paste0("10_",all_icd10$icd10)

# output time-free version of all patientid-icd10code pairs
fwrite(all_icd10,paste0("./UKB_icdcode/UKB_icd10",".txt"),sep="\t",row.names=F,col.names=T,quote=T)
##fwrite(all_icd10,paste0("./UKB_icdcode/UKB_icd10_prefix",".txt"),sep="\t",row.names=F,col.names=T,quote=T)

############################ICD 9############################
# icd9 codes with time
ukb_icd9_diagnose_time=reformatUKB(fields=c(41271, 41281),dataCoding = T)
ukb_icd9_cancer_time=reformatUKB(fields=c(40013,40005),dataCoding = T)

# icd9 diagnostic codes
icd9_diagnose_time<-ukb_icd9_diagnose_time[[1]]
colnames(icd9_diagnose_time)<-c("id","icd9","time")
# delete quotation marks in the raw data
icd9_diagnose_time<-icd9_diagnose_time[,icd9:=gsub('\\"\\"(+.+)\\"\\"',"\\1",icd9)]
icd9_diagnose_time<-icd9_diagnose_time[,icd9:=sapply(icd9,harmonizeICD9)]

# icd9 cancer codes
icd9_cancer_time<-ukb_icd9_cancer_time[[1]]
colnames(icd9_cancer_time)<-c("id","icd9","time")
icd9_cancer_time[,icd9:=sapply(icd9,harmonizeICD9)]

# combine icd9 codes from different sources
all_icd9_time<-dplyr::full_join(icd9_diagnose_time,icd9_cancer_time)

# make a NA-omitted copy for merging later
all_icd9_time_copy=na.omit(all_icd9_time)

# obtain time-free patientid-icd10code pairs data
all_icd9<-all_icd9_time[,c(1,2)]
dim(all_icd9)
# [1] 157957    2

# remove missing value
all_icd9 <- all_icd9[!(is.na(all_icd9$id)|is.na(all_icd9$icd9)),]
dim(all_icd9)
# [1] 70718     2
all_icd9 <- all_icd9[order(all_icd9$id, all_icd9$icd9),]

# ensure each patientid-icd10code pair occur but only occur once
all_icd9<-all_icd9[!duplicated(all_icd9),]
dim(all_icd9)
# [1] 70259     2

# add prefix "9_"for these icd9 codes to distinguish from icd10 codes
class(all_icd9$icd9)
#all_icd9$icd9<-paste0("09_",all_icd9$icd9)

fwrite(all_icd9,paste0("./UKB_icdcode/UKB_icd9",".txt"),sep="\t",row.names=F,col.names=T,quote=T)
##fwrite(all_icd9,paste0("./UKB_icdcode/UKB_icd9_prefix",".txt"),sep="\t",row.names=F,col.names=T,quote=T)

### merge icd9 and icd10 codes together ###

# add prefix "9_"/"10_"for these icd9/icd10 codes to distinguish from icd10/icd9 codes
#all_icd10_time_copy$icd10<-paste0("10_",all_icd10_time_copy$icd10)
#all_icd9_time_copy$icd9<-paste0("09_",all_icd9_time_copy$icd9)
# omit rows with any missing value
all_icd10_time_copy<- all_icd10_time_copy[!(is.na(all_icd10_time_copy$icd10)|is.na(all_icd10_time_copy$time)),]
all_icd9_time_copy<- all_icd9_time_copy[!(is.na(all_icd9_time_copy$icd9)|is.na(all_icd9_time_copy$time)),]
# rename columns of two datasets to the same for merge
colnames(all_icd10_time_copy)<-c("id","code","time")
colnames(all_icd9_time_copy)<-c("id","code","time")
dx_code<-rbind(all_icd10_time_copy,all_icd9_time_copy)
dim(dx_code)
# [1] 4288858       3 --- without prefix
# [1] 4288858       3 --- with prefix
## the two outcomes share the same row number, indicating no overlapping icd9/icd10 codes in UKB system,
## thus we can remain the version without the prefix for simplicity

pid <- dx_code[!duplicated(dx_code$id),"id"]
pid$Deid_ID_num <- 1:length(unique(dx_code$id))
idx <- match(dx_code$id, pid$id)
dx_code$Deid_ID_num <- pid$Deid_ID_num[idx]

dx_freq <- table(dx_code$code)
dx_unique <- names(dx_freq)

length(dx_unique)
# [1] 15204 ---with/without prefix

DxCode_num <- data.frame(DxCode=as.character(dx_unique), DxCode_num=(1:length(dx_unique)+100000000))
DxCode_num$DxCode <- as.character(DxCode_num$DxCode)

idx2 <- match(dx_code$code, DxCode_num$DxCode)
dx_code$DxCode_num <- DxCode_num$DxCode_num[idx2]

dx_code <- as.data.frame(dx_code[!is.na(idx2),c('Deid_ID_num', 'DxCode_num', 'time')])
dim(dx_code)
# [1] 4288858       3 ---with/without prefix


dx_code<-dx_code[!dx_code$time=='',]
dim(dx_code)
# [1] 4288849       3 ---with/without prefix

##dx_code<-dx_code[!duplicated(dx_code[,c("Deid_ID_num","DxCode_num")]),]
##dim(dx_code)
# [1] 4246642       3 ---with/without prefix


dx_code$DaysSinceBirth<-as.integer(difftime(dx_code$time, "1950-01-01", units = "days"))
dx_code<-dx_code[,c('Deid_ID_num','DaysSinceBirth','DxCode_num')]
dx_code<-dx_code%>%dplyr::group_by(Deid_ID_num)%>%dplyr::arrange(DaysSinceBirth,.by_group = TRUE)

write.table(DxCode_num, file='./UKB_icdcode/CodeCrosswalk/DxCode_num.txt', row.names=F)
##write.table(DxCode_num, file='./UKB_icdcode/CodeCrosswalk/DxCode_num_prefix.txt', row.names=F)
fwrite(dx_code,paste0("./UKB_icdcode/UKB_DxCode_time_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)
##fwrite(dx_code,paste0("./UKB_icdcode/UKB_DxCode_time_prefix_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)
##fwrite(dx_code,paste0("./UKB_icdcode/UKB_DxCode_time_once_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)
##tmp=fread("./UKB_icdcode/UKB_DxCode_time_once_20220622.txt",sep="\t")

## count frequency for each code (icd9 and icd10 in the same dataset)
dx_freq=as.data.frame(table(dx_code$DxCode_num))
colnames(dx_freq)=c("DxCode_num","Freq")
dx_freq$DxCode_num<-as.numeric(as.character(dx_freq$DxCode_num))
fwrite(dx_freq,paste0("./UKB_icdcode/UKB_DxCode_freq_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)

################################## Method 2: handling rare training codes ##########################
CodeAllInfo=inner_join(dx_code,DxCode_num)
CodeAllInfo=inner_join(CodeAllInfo,dx_freq)
CodeAllInfo$NewCode=''
CodeAllInfo$Source=''

threshold=5 #threshold for rare/frequent codes
FrequentCode=CodeAllInfo[CodeAllInfo$Freq>=threshold,]
FrequentCode$NewCode=FrequentCode$DxCode
FrequentCode$Source="Original"


RareCode=CodeAllInfo[!CodeAllInfo$Freq>=threshold,]
##RareCode$NewCode=ICDinDx$phecode_round[match(RareCode$DxCode,ICDinDx$ICD_dx)]

RareCode$NewCode=icd_phecode$phecode[match(RareCode$DxCode,icd_phecode$ICD_dx)]
RareCode$Source[!is.na(RareCode$NewCode)]="Phecode"
length(unique(RareCode$DxCode_num[!is.na(RareCode$NewCode)]))
# [1] 1354
# [1] 2110

RareCode$NewCode=ifelse(is.na(str_extract(RareCode$DxCode,".*\\.[0-9]")),RareCode$DxCode,RareCode$NewCode)
length(RareCode$DxCode_num[!is.na(RareCode$NewCode)])
# [1] 1365
# [1] 2893

RareCode$NewCode[is.na(RareCode$NewCode)]=ifelse(is.na(str_extract(RareCode$DxCode[is.na(RareCode$NewCode)],".*\\.[0-9]")),
                                                 RareCode$DxCode[is.na(RareCode$NewCode)],str_extract(RareCode$DxCode[is.na(RareCode$NewCode)],".*\\.[0-9]"))


RareCode$Source[(RareCode$NewCode!=RareCode$DxCode)&(RareCode$Source!="Phecode")]="Truncated"
RareCode$Source[RareCode$Source=='']="Original_Rare"
table(RareCode$Source)
sum(is.na(RareCode))


NewCodeInfo=rbind(FrequentCode,RareCode)

new_code=NewCodeInfo[,c(1,2,6,7)]

new_freq <- table(new_code$NewCode)
new_unique <- names(new_freq)

length(new_unique)
# [1] 13663
# [1] 12729
new_code<-new_code%>%dplyr::group_by(Deid_ID_num)%>%dplyr::arrange(DaysSinceBirth,.by_group = TRUE)
NewCode_num <- data.frame(NewCode=as.character(new_unique), NewCode_num=(1:length(new_unique)+100000000))
write.table(NewCode_num, file=paste0('./UKB_icdcode/CodeCrosswalk/NewCode_num_',threshold,'.txt'), row.names=F)


###############generate a crosswalk version with source#######################
NewCode_num$Source=new_code$Source[match(NewCode_num$NewCode,new_code$NewCode)]
table(NewCode_num$Source)

#      Original Original_Rare       Phecode     Truncated 
#     12316           620           524           203 
#     10903           922           642           262 
write.table(NewCode_num, file=paste0('./UKB_icdcode/CodeCrosswalk/NewCode_num_source_',threshold,'.txt'), row.names=F)
##############################################################################
idx3 <- match(new_code$NewCode, NewCode_num$NewCode)
new_code$NewCode_num <- NewCode_num$NewCode_num[idx3]

new_freq=as.data.frame(table(new_code$NewCode_num))
colnames(new_freq)=c("NewCode_num","Freq")
new_freq$NewCode_num<-as.numeric(as.character(new_freq$NewCode_num))
fwrite(new_freq,paste0("./UKB_icdcode/UKB_NewCode_freq_",today,'_',threshold,".txt"),sep="\t",row.names=F,col.names=T,quote=T)

new_code=new_code[,c(1,2,5)]
new_code<-new_code%>%dplyr::group_by(Deid_ID_num)%>%dplyr::arrange(DaysSinceBirth,.by_group = TRUE)
fwrite(new_code,paste0("./UKB_icdcode/UKB_NewCode_time_",today,'_',threshold,".txt"),sep="\t",row.names=F,col.names=T,quote=T)



