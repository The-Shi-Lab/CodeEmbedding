gradient_update_nogrp <- function(X,Y,alpha=1,convergence=1e-4){
  X = as.matrix(X)
  Y = as.matrix(Y)
  ## length normalize X and Y just in case
  X = X/apply(X,1,norm,type="2")
  Y = Y/apply(Y,1,norm,type="2")
  # Gradient update with stepsize alpha
  p = ncol(X)
  n = nrow(X)
  
  oldW <- W <- matrix(0,p,p)
  gradient = t(as.matrix(X))%*%(as.matrix(Y))
  
  error <- 10000
  while(error >=convergence){
    # Gradient update
    W <- W + alpha*gradient
    
    # perform SVD approximation
    tmp <- svd(W)
    W <- tmp$u%*%t(tmp$v)
    
    error <- sqrt(sum((oldW-W)^2))
    oldW <- W
  }## actually converges after one-step update
  return(W)
  
}

mt_method<-function(x,y,beta_hat){
  # x:p1*q, y: p2*q, beta_hat: q*q
  y_hat <- x%*%beta_hat
  y_hat <- y_hat/apply(y_hat,1,norm,type="2")
  y <- y/apply(y,1,norm,type="2")
  coss <- y_hat%*%t(y) 
  # dim(cross)=p1*p2
  
  return(coss)
}


coss_cvf = function(y_hat,y,thres, marginal, UKB_freq_norm, MGI_freq_norm, nfold){
  
  q = ncol(y_hat)
  sample.id=sample(1:nfold,size=q,replace=TRUE)
  loss = rep(NA, nfold)
  
  for(i in 1:nfold){
    # find the truncated cosine matrix
    y_hat.train = y_hat[,which(sample.id!=i),drop=FALSE]
    y.train = y[,which(sample.id!=i),drop=FALSE]
    y_hat.test = y_hat[,which(sample.id==i),drop=FALSE]
    y.test = y[,which(sample.id==i),drop=FALSE]
    
    y_hat.train.norm = y_hat.train/apply(y_hat.train, 1, norm,type="2")
    y.train.norm = y.train/apply(y.train, 1, norm,type="2")
    y_hat.test.norm = y_hat.test/apply(y_hat.test, 1, norm,type="2")
    y.test.norm = y.test/apply(y.test, 1, norm,type="2")
    
    coss = y_hat.train.norm %*% t(y.train.norm)
    if(marginal == TRUE){
      coss = refine_map_mt(map_mt = coss, source = UKB_freq_norm, target = MGI_freq_norm)
    }
    coss_trnc = coss
    coss_trnc[which(coss_trnc<thres)]=0
    
    # loss
    loss[i] = norm(y.test.norm-(t(coss_trnc) %*% y_hat.test.norm),type="F")
  }
  
  return(loss)
}

coss_cvf_leaveOneOut = function(y_hat,y,thres, marginal, UKB_freq_norm, MGI_freq_norm){
  # x:p1*q, y: p2*q, beta_hat: q*q
  # y_hat <- x%*%beta_hat
  # y_hat <- y_hat/apply(y_hat,1,norm,type="2")
  # y <- y/apply(y,1,norm,type="2")
  
  q = ncol(x)
  loss = rep(NA, q)
  for(i in 1:q){
    # find the truncated cosine matrix
    y_hat.train = y_hat[,-i]
    y.train = y[,-i]
    y_hat.test = y_hat[,i]
    y.test = y[,i]
    
    y_hat.train.norm = y_hat.train/apply(y_hat.train, 1, norm,type="2")
    y.train.norm = y.train/apply(y.train, 1, norm,type="2")
    # y_hat.test.norm = y_hat.test/apply(y_hat.test, 1, norm,type="2")
    # y.test.norm = y.test/apply(y.test, 1, norm,type="2")
    
    coss = y_hat.train.norm %*% t(y.train.norm)
    if(marginal == TRUE){
      coss = refine_map_mt(map_mt = coss, source = UKB_freq_norm, target = MGI_freq_norm)
    }
    coss_trnc = coss
    coss_trnc[which(coss_trnc<thres)]=0
    
    # loss
    loss[i] = norm(y.test-(t(coss_trnc) %*% y_hat.test),type="F")
  }
  
  return(loss)
}

coss_cv = function(y_hat,y,thres, marginal, UKB_freq_norm, MGI_freq_norm, nfold, leaveOneOut){
  loss = NULL
  for(i in 1:length(thresholds)){
    if(leaveOneOut == TRUE){
      loss.i=coss_cvf_leaveOneOut(y_hat=y_hat,y=y,thres=thresholds[i],
                               marginal=marginal, UKB_freq_norm, MGI_freq_norm)
      loss = rbind(loss,loss.i)
    }else{
      loss.i=coss_cvf(y_hat=y_hat,y=y,thres=thresholds[i],
                               marginal=marginal, UKB_freq_norm, MGI_freq_norm, nfold)
      loss = rbind(loss,loss.i)
      
    }

  }
  loss_mean = apply(loss,1,mean)
  loss_se = apply(loss,1,sd)
  loss.results = data.frame("threshold" = thresholds, "mean_loss" = loss_mean, "se" = loss_se)
  cutoff.min = thresholds[which.min(loss_mean)]
  
  return(list(loss = loss.results,
              cutoff.min = cutoff.min))
}

# find ICD codes belonging to a phecode
icd_phecode = read.csv("~/EHRmapping/mapping2/data/icd_phecodes.csv")[,-1]
icdcm_phecode = read.csv("~/EHRmapping/mapping2/data/icdcm_phecodes.csv")[,-1]

find_ICD <- function(phecode, mappingFile,poolICD, ICDversion="ICD10"){
  result = mappingFile[which(mappingFile$phecode_round == phecode),]
  if(ICDversion=="ICD9"){
    ICD_target = result$ICD_dx[!grepl("^[A-Z]",result$ICD_dx)]
  } else if(ICDversion=="ICD10"){
    ICD_target = result$ICD_dx[grepl("^[A-Z]",result$ICD_dx)]
  }
  # find ICD codes in poolICD that's in ICD_target
  selectedICD = na.omit(poolICD[poolICD %in% ICD_target])

  return(selectedICD)
}
