### Part 1 ###

# ICD to phecode crosswalk

icd9_phecode <- read.csv('/nfs/turbo/mgi-shixu/data/MGI/dx/phecode_icd9_map_unrolled.csv')
icd10_phecode <- read.csv('/nfs/turbo/mgi-shixu/data/MGI/dx/Phecode_map_v1_2_icd10_beta.csv')

colnames(icd9_phecode) <- c('ICD_dx', 'phecode')
icd9_phecode$phecode_round <- round(icd9_phecode$phecode)
icd9_phecode <- icd9_phecode[order(icd9_phecode$ICD_dx, icd9_phecode$phecode),]
icd9_phecode <- icd9_phecode[!duplicated(icd9_phecode$ICD_dx),]

icd10_phecode <- icd10_phecode[,c('ICD10','PHECODE')]
colnames(icd10_phecode) <- c('ICD_dx', 'phecode')
icd10_phecode$phecode_round <- round(icd10_phecode$phecode)
icd10_phecode <- icd10_phecode[order(icd10_phecode$ICD_dx, icd10_phecode$phecode),]
icd10_phecode <- icd10_phecode[!duplicated(icd10_phecode$ICD_dx),]

if(Version=='ICD9'){icd_phecode <- icd9_phecode}
if(Version=='ICD10'){icd_phecode <- icd10_phecode}
if(Version=='ICD9ICD10'){icd_phecode <- rbind(icd9_phecode, icd10_phecode)}

# in SurgPcpCohort #
#ICD9: 7708
#ICD10: 3424
#ICD9ICD10: 11132




### Part 2 ###

# code to calculate cosine similarty
# CodeEmbedVec_tmp: code embedding matrix for 500 vectors (formatted to match phecode outcome)

for(j in c(10,30,50,100,150,200,250,300,350,400,450,500)){
  tmp <- CodeEmbedVec_tmp[,1:j]
  tmp2 <- sqrt(apply(tmp^2,1,sum))
  tmp_norm <- tmp/tmp2
  cosine_matrix <- tmp_norm%*%t(tmp_norm)
  cosine_vec <- cosine_matrix[upper.tri(cosine_matrix)]

  print(j)
}


