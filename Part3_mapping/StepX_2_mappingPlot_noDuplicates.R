library(readxl)
library(tidyr)
library(ggplot2)
library(dplyr)

# orgCos = read_excel("/Users/jiacong/Google Drive/Umich/research/EHRmapping/mapping/results/cos_target_401_test/cos_target_401_totalObs_noDuplicates/Coss_target_ICD10_SR_pmi_vec500_refined.xlsx", col_names=FALSE)
# truncatedCos = read_excel("/Users/jiacong/Google Drive/Umich/research/EHRmapping/mapping/results/cos_target_401_test/cos_target_401_totalObs_noDuplicates/Coss_target_ICD10_SR_pmi_vec500_refined_truncated.xlsx", col_names=FALSE)
# title = "ICD10 mapping from MGI to UKB using the updated cosine (without duplicates)"
# file_name ="/Users/jiacong/Google Drive/Umich/research/EHRmapping/mapping/results/cos_target_401_test/cos_target_401_totalObs_noDuplicates/ICD10_updated.png"

# MGI_embed_duplicates = FALSE
phecode = 401
outpath = "/net/snowwhite/home/jiacong/EHRmapping/mapping2/results/"
cosineType = "orgCosine" # refinedCosine
cosineType = "refinedCosine"

# orgCos = read.csv(paste0("~/EHRmapping/mapping2/results/phecode_401_",cosineType,"_withDuplicates.csv"),header = TRUE, row.names = 1)
# truncatedCos = read.csv(paste0("~/EHRmapping/mapping2/results/phecode_401_",cosineType,"_withDuplicates_thres.csv"),header = TRUE, row.names = 1)
# title = "ICD10 mapping from MGI to UKB using the refined cosine (with duplicates)"
# file_name ="~/EHRmapping/mapping2/results/phecode_401_withDuplicates_refined.png"

orgCos = read.csv(paste0("~/EHRmapping/mapping2/results/phecode_401_",cosineType,"_noDuplicates.csv"),header = TRUE, row.names = 1)
truncatedCos = read.csv(paste0("~/EHRmapping/mapping2/results/phecode_401_",cosineType,"_noDuplicates_thres.csv"),header = TRUE, row.names = 1)
title = "ICD10 mapping from MGI to UKB using the refined cosine (without duplicates)"
file_name ="~/EHRmapping/mapping2/results/phecode_401_noDuplicates_refined.png"


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

cbind(MGI_codes,MGI_freq)
cbind(UKB_codes,UKB_freq)


