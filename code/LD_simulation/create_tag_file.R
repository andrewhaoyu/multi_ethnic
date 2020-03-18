#create tag file for LD simulation
#tag SNPs are Biallelic_SNP with MAF 1% in EUR, AFR, AMR, EAS and SAS
#tag files are created seperately for each chr
library(data.table)
library(dplyr)
eth <- c("AFR","AMR","EAS","EUR","SAS")
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  TYPE = leg$TYPE
  MAF <- leg[,c(6:10)]
  for(k in 1:length(eth)){
    idx <- which(MAF[,k]>=0.01&
                   MAF[,k]<=0.99&
                   TYPE=="Biallelic_SNP")
    leg_chr = leg$position[idx] 
    
    
    write.table(leg_chr,file = paste0("/data/zhangh24/KG.impute2/tag/",eth[k],"_chr",i,".tag"),row.names = F,col.names = F,quote=F)
    
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
for(i in 1:4){
  for(j in 1:22){
    
    tag <- tag_all[[j]]
    gen <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",j,"_",k,".controls.gen"),header=F))
    colnames(tag) <- "position"
    colnames(gen)[3] <- "position"
    tag.gen <- left_join(tag,gen,by="position")
    #reorder the column back to the original order
    tag.gen <- tag.gen[,c(2,3,1,4,5)]
    write.table(tag.gen,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag.info.txt"),quote = F,col.names = F,row.names = F)
  }
  
  
  
}
  







