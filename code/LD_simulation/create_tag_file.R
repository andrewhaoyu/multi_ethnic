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


eth <- c("EUR","AFR","AMR","EAS")
#read the genotype data
#filter the data in the tag file
#i is the ethnic group
#j is the chr
#k is the replicates for each job
#in gen file
#create tag info file first
library(dplyr)
for(i1 in 1:1){
  for(j in 1:22){
    #since the SNP code is the same, just need the first one
    k <- 1
    tag <- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i1],"_chr",j,".tag"),header=F)
    
    gen <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",j,"_",k,".cases.gen"),header=F))
    colnames(tag) <- "position"
    colnames(gen)[3] <- "position"
    tag.gen <- left_join(tag,gen,by="position")
    #reorder the column back to the original order
    tag.info <- tag.gen[,c(2,3,1,4,5)]
    write.table(tag.info,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",j,".tag.info.txt"),quote = F,col.names = F,row.names = F)
  }
  
  
  
}







#count the number of SNPs in each ethnic group
eth <- c("EUR","AFR","AMR","EAS")
SNP.number <- matrix(0,22,4)
for(i1 in 1:4){
  for(i2 in 1:22){
    file <- paste0(" /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",i2,".tag.info.txt")  
    code <- paste0("wc -l",file)
    SNP.number[i2,i1] <- as.numeric(gsub(file,"",system(code,intern=T)))
  }
}




#generate all SNPs information list
library(data.table)
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS")
snp.infor.list <- list()
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  TYPE = leg$TYPE
  MAF <- leg[,c(9,6,7,8)]
  MAF.Bi <- MAF
  for(k in 1:4){
  MAF.Bi[,k] <- ifelse(MAF[,k]>=0.01&
                     MAF[,k]<=0.99,1,0)
  }
    idx <- which((MAF.Bi[,1]==1|
                    MAF.Bi[,2]==1|
                    MAF.Bi[,3]==1|
                    MAF.Bi[,4]==1)&
              (TYPE=="Biallelic_SNP"))
    CHR <- rep(i,length(idx))
    snp.infor.temp <- cbind(leg[idx,],CHR)
    snp.infor.list[[i]] <- snp.infor.temp
    
}

snp.infor <- rbindlist( snp.infor.list)
save(snp.infor,file =  "/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
