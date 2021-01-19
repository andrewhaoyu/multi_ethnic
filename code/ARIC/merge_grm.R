for(i in 1:2){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
  temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[1],"/",eth[i],"/")
  data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
  out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
  
  multi_grm = rep("c",22)
  for(j in 1:22){
    multi_grm[j] = paste0(temp.dir,"chr.qc",j)
  }
  write.table(multi_grm,file =paste0(temp.dir,"multi_grm.txt"),col.names = F,row.names = F,quote=F)
  res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/gcta/gcta64 --mgrm ",temp.dir,"multi_grm.txt --make-grm --out ",temp.dir,"merged_chr_qc"))  
}
