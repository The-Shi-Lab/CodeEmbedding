library(ggplot2)
library(reshape2)
##### import pvalues for code-wise comparison #####
load("~/EHRmapping/comparison/results/codeWise_step1_1_data_ICD9_step2.RData")
rm(list=setdiff(ls(),"pvalues.codeWise.df"))

#### import pvalues for block-wise comparison ####
load("~/EHRmapping/comparison/results/blockWise_step1_2_ICD9_analysis.RData")

pvalues.burden = pvalues.burden.df
pvalues.skat = pvalues.skat.df

##### 
code_block = block_code_crosswalk
all_blocks = unique(code_block$block)
all_codes = unique(block_code_crosswalk$code)

block_new = data.frame("block" = all_blocks,"block_new" = 1:length(all_blocks))
code_block$block_new = block_new$block_new[match(code_block$block, block_new$block)]
pvalues_results = data.frame("codes" = all_codes, 
                             "pvalues_codes"=pvalues.codeWise.df$pvalue[match(all_codes,pvalues.codeWise.df$code)],
                             "block"= code_block$block[match(all_codes,code_block$code)],
                             "block_new"= code_block$block_new[match(all_codes,code_block$code)])
pvalues_results2 = pvalues_results[order(pvalues_results$block_new),]
pvalues_results2 = na.omit(pvalues_results2)

pvalues.burden.df = data.frame('block'=all_blocks,pvalues=pvalues.burden$pvalue[match(all_blocks,pvalues.burden$block)])
pvalues.burden.df$block_new = block_new$block_new[match(pvalues.burden.df$block, block_new$block)]
pvalues.skat.df = data.frame('block'=all_blocks,pvalues=pvalues.skat$pvalue[match(all_blocks,pvalues.skat$block)])
pvalues.skat.df$block_new = block_new$block_new[match(pvalues.skat.df$block, block_new$block)]

pvalues.block = pvalues.burden.df
pvalues.block$SKAT  = pvalues.skat.df$pvalues[match(pvalues.block$block,pvalues.skat.df$block)]
colnames(pvalues.block)  = c("block","Burden","block_new","SKAT")


##### make the table #####
load("~/EHRmapping/comparison/results/blockWise_step1_1_ICD9_prepareData.RData")

p.threshold = 0.05/length(all_blocks)
## block-burden
sig.burden.ind = pvalues.burden.df$block[which(pvalues.burden.df$pvalue<=p.threshold)]
length(sig.burden.ind) 

## block-skat
sig.skat.ind = pvalues.skat.df$block[which(pvalues.skat.df$pvalue<=p.threshold)]
length(sig.skat.ind) 

block_indx = setdiff(sig.skat.ind, sig.burden.ind)
length(block_indx)

df=NULL
UKB_id = unique(UKB_ICD9_joint$id)
MGI_id = unique(MGI_ICD9_joint$id)
for(i in 1:length(block_indx)){
  print(i)
  block=block_indx[i]
  p.burden = pvalues.burden.df$pvalues[which(pvalues.burden.df$block==block)]
  p.skat = pvalues.skat.df$pvalues[which(pvalues.skat.df$block==block)]
  codes = block_code_crosswalk$code[which(block_code_crosswalk$block==block)]
  
  MGI_ICD9_wide = UKB_ICD9_wide = NULL
  for(j in 1:length(codes)){
    code = codes[j]
    
    selectedID = unique(MGI_ICD9_joint$id[which(MGI_ICD9_joint$icd9 == code)])
    temp = (MGI_id %in% selectedID)*1
    MGI_ICD9_wide = cbind(MGI_ICD9_wide, temp)
    
    selectedID = unique(UKB_ICD9_joint$id[which(UKB_ICD9_joint$icd9 == code)])
    temp = (UKB_id %in% selectedID)*1
    UKB_ICD9_wide = cbind(UKB_ICD9_wide, temp)
    
  }
  UKB_data = as.matrix(UKB_ICD9_wide)
  MGI_data = as.matrix(MGI_ICD9_wide)
  
  df_temp = data.frame("UKB"=apply(UKB_data,2,sum), "MGI"= apply(MGI_data,2,sum))
  rownames(df_temp) = codes
  total = data.frame("subtotal"=apply(df_temp,2,sum))
  df_temp = rbind(df_temp,t(total))
  UKB_norm = df_temp$UKB/nrow(UKB_data)
  MGI_norm = df_temp$MGI/nrow(MGI_data)
  df_temp = cbind(df_temp, UKB_norm, MGI_norm)
  df_temp = cbind(block,p.skat,p.burden,df_temp)
  
  df= rbind(df,df_temp)
}
head(df)

write.csv(df,"~/EHRmapping/comparison/results/step2_2_ICD9_skat_vs_burden_pBonferroni.csv")
# write.csv(df,"~/EHRmapping/comparison/results/step2_2_ICD9_skat_vs_burden_p0.csv")



##### draw the plot ######

# make 0 = 0+epsilon
pvalue.codeWise.min = min(pvalues_results2$pvalues_codes[-which(pvalues_results2$pvalues_codes==0)])
pvalue.burden.min = min(pvalues.block$Burden[-which(pvalues.block$Burden==0)])
pvalue.skat.min = min(pvalues.block$SKAT[-which(pvalues.block$SKAT==0)])
pvalue.min = 10^(log10(min(c(pvalue.codeWise.min,pvalue.burden.min, pvalue.skat.min)))-0.5)
pvalues_results2$pvalues_codes[which(pvalues_results2$pvalues_codes==0)] = pvalue.min
pvalues.block$Burden[which(pvalues.block$Burden==0)] = pvalue.min
pvalues.block$SKAT[which(pvalues.block$SKAT==0)] = pvalue.min

pvalues.block.long = melt(pvalues.block[,-1],id.vars="block_new")
head(pvalues.block.long)
colnames(pvalues.block.long) = c("block_new","Method","pvalues")

block_thres = 0.05/length(all_blocks)


(pic = ggplot(data=pvalues_results2, aes(x=block_new,y=-log10(pvalues_codes),color=block_new))+
    geom_point()+
    scale_color_gradient2(low = "yellow", mid = "darkblue", high = "red", midpoint = mean(pvalues_results2$block_new), guide="none")+
    geom_point(data=pvalues.block.long,aes(x=block_new,y=-log10(pvalues),shape=Method),size=1.5,color="red")+
    scale_shape_manual(values=c(2,4))+
    geom_hline(yintercept = -log10(block_thres))+
    xlab("Block")+
    ylab("-log10(pvalue)")+
    ggtitle("ICD-9 codes")+
    theme_bw()+
    theme(axis.text.x=element_blank()) )

ggsave(pic,file="~/EHRmapping/comparison/results/step2_1_ICD9_compare_MhPlot.pdf",device = "pdf",width=16,height = 8)




