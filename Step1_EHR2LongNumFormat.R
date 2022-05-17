

### preprocess dataset into co-occurrence format ###



## 1)
## create ID crosswalk for patient ID, dx, px, and lab ##


# Pt ID #

library(data.table)
demo <- fread('/nfs/turbo/mgi-shixu/data/MGI/raw/HPI_4605_Demographics.txt')
demo <- as.data.frame(demo)
demo$Deid_ID_num <- 1:nrow(demo)

# write.table(demo, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/PtID_num.txt', row.names=F)


# diagnosis #

dx_code <- fread('/nfs/turbo/mgi-shixu/data/MGI/raw/HPI_4605_Diagnosis.txt', select=3)
dx_freq <- table(dx_code$DiagnosisCode)
idx <- dx_freq>10
dx_unique <- names(dx_freq[idx])

length(dx_unique)
# [1] 27014

DxCode_num <- data.frame(DxCode=as.character(dx_unique), DxCode_num=(1:length(dx_unique)+100000000))
DxCode_num$DxCode <- as.character(DxCode_num$DxCode)

# write.table(DxCode_num, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/DxCode_num.txt', row.names=F)


# procedure #

px_code <- fread('/nfs/turbo/mgi-shixu/data/MGI/raw/HPI_4605_Procedures.txt', select=3)
px_freq <- table(px_code$ProcedureCode)
idx <- px_freq>10
px_unique <- names(px_freq[idx])

length(px_unique)
# [1] 8719

PxCode_num <- data.frame(PxCode=as.character(px_unique), PxCode_num=(1:length(px_unique)+200000000))
PxCode_num$PxCode <- as.character(PxCode_num$PxCode)

# write.table(PxCode_num, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/PxCode_num.txt', row.names=F)


# lab #

lab_code <- fread('/nfs/turbo/mgi-shixu/data/MGI/lab_abn/AllAbnLabResults.txt.gz', select=6)
# all abn lab results were filter by min count of 10 previously
lab_unique <- unique(lab_code$ResultCode)

length(lab_unique)
# [1] 603

LabCode_num <- data.frame(LabCode=as.character(lab_unique), LabCode_num=(1:length(lab_unique)+300000000))
LabCode_num$LabCode <- as.character(LabCode_num$LabCode)

# write.table(LabCode_num, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/LabCode_num.txt', row.names=F)






## 2)
## make dx, px, and lab into all numeric variables
## remove records with missing DaysSinceBirth

library(data.table)
PtID_Crosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/PtID_num.txt', header=T)


# diagnosis
dx <- fread('/nfs/turbo/mgi-shixu/data/MGI/raw/HPI_4605_Diagnosis.txt', select=c(1,2,3))
colnames(dx) <- c('Deid_ID', 'DaysSinceBirth', 'Code')
dx_crosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/DxCode_num.txt', header=T)
dim(dx)
# [1] 118347476         3

idx <- match(dx$Deid_ID, PtID_Crosswalk$Deid_ID)
dx$Deid_ID_num <- PtID_Crosswalk$Deid_ID_num[idx]

idx2 <- match(dx$Code, dx_crosswalk$DxCode)
dx$Code_num <- dx_crosswalk$DxCode_num[idx2]

dx_sort <- as.data.frame(dx[!is.na(idx2),c('Deid_ID_num', 'DaysSinceBirth', 'Code_num')])
dim(dx_sort)
# [1] 118295661         3
dx_sort <- dx_sort[!is.na(dx_sort$DaysSinceBirth),]
dim(dx_sort)
# [1] 118072880         3
dx_sort <- dx_sort[order(dx_sort$Deid_ID_num, dx_sort$DaysSinceBirth),]

write.csv(dx_sort, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/LongFormat_Num/Dx_num.csv', row.names=F)
# save(dx_sort, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/LongFormat_Num/Dx_num.RData')



# procedure
px <- fread('/nfs/turbo/mgi-shixu/data/MGI/raw/HPI_4605_Procedures.txt', select=c(1,2,3))
colnames(px) <- c('Deid_ID', 'DaysSinceBirth', 'Code')
px_crosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/PxCode_num.txt', header=T)
dim(px)
# [1] 10063101        3

idx <- match(px$Deid_ID, PtID_Crosswalk$Deid_ID)
px$Deid_ID_num <- PtID_Crosswalk$Deid_ID_num[idx]

idx2 <- match(px$Code, px_crosswalk$PxCode)
px$Code_num <- px_crosswalk$PxCode_num[idx2]

px_sort <- as.data.frame(px[!is.na(idx2), c('Deid_ID_num','DaysSinceBirth','Code_num')])
dim(px_sort)
# [1] 10026348        3
px_sort <- px_sort[!is.na(px_sort$DaysSinceBirth),]
dim(px_sort)
# [1] 10026347        3
px_sort <- px_sort[order(px_sort$Deid_ID_num, px_sort$DaysSinceBirth),]

write.csv(px_sort, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/LongFormat_Num/Px_num.csv', row.names=F)
# save(px_sort, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/LongFormat_Num/Px_num.RData')



# lab
lab <- fread('/nfs/turbo/mgi-shixu/data/MGI/lab_abn/AllAbnLabResults.txt.gz', select=c(1,2,6))
colnames(lab) <- c('Deid_ID', 'DaysSinceBirth', 'Code')
lab_crosswalk <- read.table('/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/CodeCrosswalk/LabCode_num.txt', header=T)
dim(lab)
# [1] 9370966        3

idx <- match(lab$Deid_ID, PtID_Crosswalk$Deid_ID)
lab$Deid_ID_num <- PtID_Crosswalk$Deid_ID_num[idx]

idx2 <- match(lab$Code, lab_crosswalk$LabCode)
lab$Code_num <- lab_crosswalk$LabCode_num[idx2]

lab_sort <- as.data.frame(lab[!is.na(idx2), c('Deid_ID_num','DaysSinceBirth','Code_num')])
dim(lab_sort)
# [1] 9370966        3
lab_sort <- lab_sort[!is.na(lab_sort$DaysSinceBirth),]
dim(lab_sort)
# [1] 9370966        3
lab_sort <- lab_sort[order(lab_sort$Deid_ID_num, lab_sort$DaysSinceBirth),]

write.csv(lab_sort, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/LongFormat_Num/Lab_num.csv', row.names=F)
# save(lab_sort, file='/nfs/turbo/mgi-shixu/project/CodeEmbedding/data/LongFormat_Num/Lab_num.RData')

