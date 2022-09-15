### clean up the results ###
library(ggplot2)
DataSources <- 'New'
#DataSources <- 'Dx'
Versions <- c('ICD9ICD10')
n_vec <- c(10,30,50,100,150,200,250,300,350,400,450,500)
#Windows <- c(1, 7, 10, 14, 20, 30, 40, 50, 60)
Windows <- c(0, 1, 6, 9, 13, 19, 29, 39, 49, 59)
threshold <-2

file_path <- paste0('/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/NewAUCsummary/',threshold,'/')
#file_path <- '/net/junglebook/home/shuboz/createUKBphenome/UKB_icdcode/Evaluation/DxAUCsummary/'

for(DataSource in DataSources){
  for(Version in Versions){
    out <- data.frame(row.names=paste0('window_', Windows, 'days'))
    out2 <- data.frame(row.names=paste0('window_', Windows, 'days'))
    for(i in n_vec){
      file_name <- paste0('AUC_', DataSource, i, 'vec_', Version, '.csv')
      tmp <- read.csv(paste0(file_path, file_name))
      tmp2 <- round(tmp, 4)
      auc <- formatC(tmp2[,1], 4, format='f')
      lcl <- formatC(tmp2[,2], 4, format='f')
      ucl <- formatC(tmp2[,3], 4, format='f')
      out_tmp <- paste0(auc, ' (', lcl, '-', ucl, ')')
      out <- cbind(out, out_tmp)
      out2<-cbind(out2,auc)
    }
    colnames(out) <- paste0('emb_', n_vec, 'vec')
    colnames(out2) <- paste0('emb_', n_vec, 'vec')
    write.csv(out, file=paste0(file_path, DataSource, '_' , Version, '_summary.csv'), row.names=F)
    write.csv(out2, file=paste0(file_path, DataSource, '_' , Version, '_summary2.csv'), row.names=F)
  }
}

out2<-as.matrix(apply(out2, 2, as.numeric))
rownames(out2)<-paste0("W",Windows)
out2_melt<-melt(out2)
ggp <- ggplot(out2_melt, aes(Var1, Var2,fill=value)) + geom_tile()+ scale_fill_distiller(palette = "RdPu")
ggp                   
