library(gghighlight)
library(ggplot2)
library(reshape)
library(ggpubr)
library(gridExtra)
library(grid)

DataSources <- c('New')
Versions <- c('ICD9ICD10')
n_vec <- c(10,30,50,100,150,200,250,300,350,400,450,500)
Windows <- c(0, 1, 6, 9, 13, 19, 29, 39, 49, 59)+1

file_path <- '/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/Evaluation/NewAUCsummary/'

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

out2_<-as.matrix(apply(out2, 2, as.numeric))
rownames(out2_)<-as.character(paste0(Windows))
out2_
out2_melt<-melt(out2_)
colnames(out2_melt) = c('Windows','Embed_vec','AUC')
head(out2_melt)
out2_melt$Embed_vec<-gsub(".*_","",out2_melt$Embed_vec)
out2_melt$Embed_vec<-as.character(gsub("vec","",out2_melt$Embed_vec))
out2_melt$Embed_vec<-factor(out2_melt$Embed_vec,levels=c("10","30","50","100","150","200","250","300","350","400","450","500"))
out2_melt$Windows<-factor(out2_melt$Windows,levels=c("1","2","7","10","14","20","30","40","50","60"))

loc = which(out2_ == max(out2_), arr.ind = TRUE)
(pic = ggplot(out2_melt, aes(x=Windows, y=Embed_vec,fill=AUC)) + 
  geom_tile()+ 
  xlab('Time window in days')+
  ylab('Length of embedding vectors')+
  ggtitle('AUC of embeddings for main diagnosis in MGI (with repeated measurements)')+
  scale_fill_distiller(palette = "RdPu",limits=c(0.5,0.9),direction="reverse") +
  geom_point(data = out2_melt[which.max(out2_melt$AUC), ], color="red", 
             size=3) +
  geom_text(data = out2_melt[which.max(out2_melt$AUC), ], aes(x=loc[1]+0.5,y=loc[2]+0.5,label=max(out2_melt$AUC)))+
  theme_bw()+
  theme(text = element_text(size = 12),
        legend.position = "bottom") +    
  labs(tag="B") )
pic_MGI = pic

ggsave('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/MGI_icdcode_thres10_AUC.pdf', device = "pdf", width = 8, height = 8)

pic1=ggarrange(pic_UKB,pic_MGI,nrow=1)
pic1
# combine pic and pic_UKB into one figure
ggsave('/net/snowwhite/home/jiacong/EHRmapping/MGI_embedding/MGI_thres10_mainDiagnosis/AUC_combined.pdf', device = "pdf", width = 14, height = 7)


