library(dplyr)

# Signifiant unique blocks

load("~/EHRmapping/comparison/results/codeWise_step1_1_data_ICD10_step2.RData")
rm(list=setdiff(ls(),"pvalues.codeWise.df"))

#### import pvalues for block-wise comparison ####
load("~/EHRmapping/comparison/results/blockWise_step1_4_ICD10_analysiss.RData")

pvalues.burden = pvalues.burden.df
pvalues.skat = pvalues.skat.df

n.block.total = length(unique(pvalues.burden$block))


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

######### significant code ##########
nCodes = length(unique(pvalues_results2$codes))
n_p0 = sum(pvalues_results2$pvalues_codes<=0)
prop_p0 = mean(pvalues_results2$pvalues_codes<=0)

threshold_code_BH = 0.05/nCodes
n_pB = sum(pvalues_results2$pvalues_codes<=threshold_code_BH)
prop_pB = mean(pvalues_results2$pvalues_codes<=threshold_code_BH)

codeWise_results = rbind(paste0(n_p0," (", round(prop_p0,2), ") "),
                         paste0(n_pB," (", round(prop_pB,2), ") ") )
codeWise_results = cbind(c(nCodes,nCodes),codeWise_results)
rownames(codeWise_results) = c("p0","pB")
codeWise_results
write.csv(codeWise_results, "~/EHRmapping/comparison/results/step3_2_ICD10_numberOfsigCodes.csv")


########### significant blocks ######
# p=0
p.threshold = 0
length(pvalues.burden.df$pvalues)

## block-burden
sig.burden.ind = pvalues.burden.df$block[which(pvalues.burden.df$pvalue<=p.threshold)]
length(sig.burden.ind) 

## block-skat
sig.skat.ind = pvalues.skat.df$block[which(pvalues.skat.df$pvalue<=p.threshold)]
length(sig.skat.ind) 

sig.unique.skat = setdiff(sig.skat.ind, sig.burden.ind)
length(sig.unique.skat)

sig.unique.burden = setdiff(sig.burden.ind, sig.skat.ind)
length(sig.unique.burden)
sig.summary.p0 = data.frame("Total"=n.block.total,"Burden"=length(sig.burden.ind),
                            "SKAT"=length(sig.skat.ind),"Burden\\SKAT" = length(sig.unique.burden),
                            "SKAT\\Burden" = length(sig.unique.skat))

# p=BH
p.threshold = 0.05/length(all_blocks)
# p.threshold = 0

## block-burden
sig.burden.ind = pvalues.burden.df$block[which(pvalues.burden.df$pvalue<=p.threshold)]
length(sig.burden.ind) 

## block-skat
sig.skat.ind = pvalues.skat.df$block[which(pvalues.skat.df$pvalue<=p.threshold)]
length(sig.skat.ind)

sig.unique.skat = setdiff(sig.skat.ind, sig.burden.ind)
length(sig.unique.skat)

sig.unique.burden = setdiff(sig.burden.ind, sig.skat.ind)
length(sig.unique.burden)
sig.summary.pBH = data.frame("Total"=n.block.total,"Burden"=length(sig.burden.ind),
                             "SKAT"=length(sig.skat.ind),"Burden-SKAT" = length(sig.unique.burden),
                             "SKAT-Burden" = length(sig.unique.skat))

sig.summary = rbind(sig.summary.p0, sig.summary.pBH)
rownames(sig.summary) = c("p=0","p<p.Bonferroni")


sig.summary2 = sig.summary %>%
  mutate(Burden1 = paste(Burden," (",round(Burden/Total,2),")",sep="")) %>%
  mutate(SKAT1 = paste(SKAT," (",round(SKAT/Total,2),")",sep="")) %>%
  mutate(Burden.SKAT1 = paste(Burden.SKAT," (",round(Burden.SKAT/Total,2),")",sep="")) %>%
  mutate(SKAT.Burden1 = paste(SKAT.Burden," (",round(SKAT.Burden/Total,2),")",sep="")) %>%
  select(c(Total,Burden1, SKAT1, Burden.SKAT1, SKAT.Burden1))

sig.summary2

write.csv(sig.summary2, "~/EHRmapping/comparison/results/step3_2_ICD10_numberOfsigBlocks.csv")
