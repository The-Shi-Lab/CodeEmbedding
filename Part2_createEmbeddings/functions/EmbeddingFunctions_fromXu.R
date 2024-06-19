library(gmp)
library(Matrix)
library(irlba)
library(readr)
library(data.table)
library(dplyr)
options(stringsAsFactors = FALSE)

construct_pmi <- function(coccur,singletons,my.smooth=0.75){
  names(coccur) = c("code1","code2","joint_count")
  ind <- which(coccur$code1!=coccur$code2 &
                 coccur$code1%in%singletons$marg_word$code & 
                 coccur$code2%in%singletons$marg_context$code)
  coccur = coccur[ind,]
  
  coccur$joint_count = as.numeric(coccur$joint_count)
  singletons$marg_word$marg_count = as.numeric(singletons$marg_word$marg_count)
  singletons$marg_context$marg_count = as.numeric(singletons$marg_context$marg_count)
  
  pmi_df <- coccur %>%
    inner_join(singletons$marg_word,by=c("code1" = "code")) %>%
    dplyr::rename(W=marg_count) %>%
    inner_join(singletons$marg_context,by=c("code2" = "code")) %>%
    dplyr::rename(C=marg_count) %>%
    mutate(  PMI = joint_count/(W * ((C/singletons$D)^my.smooth))  ) %>%
    mutate(PMI=log(PMI)) 
  # %>%select(code1,code2,PMI)
  return(pmi_df)
}

construct_sppmi <- function(pmi,k=10) {
  sppmi_df <- pmi %>%
    mutate(SPPMI = pmax(PMI - log(k),0))%>%  
    select(code1,code2,SPPMI)
  
  # sppmi_df <- pmi %>%
  #   mutate(SPPMI = PMI) 
  
  all_words <- unique(c(sppmi_df$code1,sppmi_df$code2))
  word_2_index <- 1:length(all_words)
  names(word_2_index) <- all_words
  
  i <- as.numeric(word_2_index[as.character(sppmi_df$code1)])
  j <- as.numeric(word_2_index[as.character(sppmi_df$code2)])
  x <- as.numeric(sppmi_df$SPPMI)
  
  ## Remove 0s ##
  non_zero <- which(x != 0)
  i <- i[non_zero]
  j <- j[non_zero]
  x <- x[non_zero]
  if(max(i)<length(all_words)|max(j)<length(all_words)){
    i=c(i,length(all_words))
    j=c(j,length(all_words))
    x=c(x,0)
  }
  
  ism <- c(i,j)
  jsm <- c(j,i)
  xsm <- c(x,x)
  sppmi <- sparseMatrix(i=ism,j=jsm,x=xsm)
  rownames(sppmi) <- all_words
  colnames(sppmi) <- all_words
  return(sppmi)
}

construct_sppmi_pmi <- function(pmi,k=10) {
  # sppmi_df <- pmi %>%
  #   mutate(SPPMI = pmax(PMI - log(k),0)) 
  # %>%  select(code1,code2,SPPMI)
  
  sppmi_df <- pmi %>%
    mutate(SPPMI = PMI) 
  
  all_words <- unique(c(sppmi_df$code1,sppmi_df$code2))
  word_2_index <- 1:length(all_words)
  names(word_2_index) <- all_words
  
  i <- as.numeric(word_2_index[as.character(sppmi_df$code1)])
  j <- as.numeric(word_2_index[as.character(sppmi_df$code2)])
  x <- as.numeric(sppmi_df$SPPMI)
  
  ## Remove 0s ##
  non_zero <- which(x != 0)
  i <- i[non_zero]
  j <- j[non_zero]
  x <- x[non_zero]
  if(max(i)<length(all_words)|max(j)<length(all_words)){
    i=c(i,length(all_words))
    j=c(j,length(all_words))
    x=c(x,0)
  }
  
  ism <- c(i,j)
  jsm <- c(j,i)
  xsm <- c(x,x)
  sppmi <- sparseMatrix(i=ism,j=jsm,x=xsm)
  rownames(sppmi) <- all_words
  colnames(sppmi) <- all_words
  return(sppmi)
}

factor_sppmi <- function(sppmi,dim_size=100,iters=25,remove_empty=TRUE,use_sum=F) {
  fit <- irlba(sppmi,nv=dim_size,maxit=iters,verbose=TRUE)
  W <- fit$u %*% diag(sqrt(fit$d))
  vecs <- W 
  if(use_sum) {
    C <- fit$v %*% diag(sqrt(fit$d))  
    vecs <- vecs + C
  }
  rownames(vecs) <- rownames(sppmi)
  if(remove_empty) {
    ## Remove empty word vectors ##
    vecs <- vecs[which(rowSums(abs(vecs)) != 0),]
  }
  return(list(vecs=vecs,fit=fit))
}

getSingletonTb <- function(coccur){
  marg_word = coccur %>% group_by(code1) %>% summarise(marg=sum(count))
  marg_context = coccur %>% group_by(code2) %>% summarise(marg=sum(count))
  names(marg_word) = c("code","marg_count")
  names(marg_context) = c("code","marg_count")
  D = sum(as.numeric(coccur$count))
  return(list(marg_word=marg_word,marg_context=marg_context,D=D))
}
