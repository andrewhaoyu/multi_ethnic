#create tag file for LD simulation
#tag SNPs are Biallelic_SNP with MAF 1% in EUR, AFR, AMR, EAS and SAS
#tag files are created seperately for each chr
library(data.table)
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
for(i in 1:22){
  print(i)
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  TYPE = leg$TYPE
  MAF <- leg[,c(9,6,7,8,10)]

  for(k in 1:length(eth)){
  if(k==1){
    idx <- which(MAF[,k]>=0.01&
                   MAF[,k]<=0.99&
                   TYPE=="Biallelic_SNP")
    leg_chr = leg$position[idx] 
    
    
    write.table(leg_chr,file = paste0("/data/zhangh24/KG.impute2/tag/",eth[k],"_chr",i,".tag"),row.names = F,col.names = F,quote=F)
    
  }else{
    idx <- which(MAF[,k]>=0.005&
                   MAF[,k]<=0.995&
                   TYPE=="Biallelic_SNP")
    leg_chr = leg$position[idx] 
    
    
    write.table(leg_chr,file = paste0("/data/zhangh24/KG.impute2/tag/",eth[k],"_chr",i,".tag"),row.names = F,col.names = F,quote=F)
    
  }
    
  }
  
  
}


#generate all SNPs information list
library(data.table)
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
snp.infor.list <- list()
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  TYPE = leg$TYPE
  MAF <- leg[,c(9,6,7,8,10)]
  MAF.Bi <- MAF
  for(k in 1:length(eth)){
  MAF.Bi[,k] <- ifelse(MAF[,k]>=0.01&
                     MAF[,k]<=0.99,1,0)
  }
    idx <- which((MAF.Bi[,1]==1|
                    MAF.Bi[,2]==1|
                    MAF.Bi[,3]==1|
                    MAF.Bi[,4]==1|
                    MAF.Bi[,5]==1)&
              (TYPE=="Biallelic_SNP"))
    CHR <- rep(i,length(idx))
    snp.infor.temp <- cbind(leg[idx,],CHR)
    snp.infor.list[[i]] <- snp.infor.temp
    
}

snp.infor <- rbindlist(snp.infor.list)
save(snp.infor,file =  "/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
