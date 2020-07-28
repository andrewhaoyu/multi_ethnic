#goal extract causal SNPs

#load causal SNPs
setwd("/data/zhangh24/multi_ethnic")
load("./result/LD_simulation/cau.snp.infor.rdata")
#write the target SNPs id into the file
snp.id <- cau.snp.infor$id
write.table(snp.id,file  ="./result/LD_simulation/causal_snp_id.txt",quote = F,row.names = F,col.names = F )

#write the target SNPs by chr and position
for(j in 1:22){
  jdx <- which(cau.snp.infor$CHR==j)
  snp.pos <-cau.snp.infor[jdx,]$position
  write.table(snp.pos,file  =paste0("./result/LD_simulation/causal_snp_pos_chr_",j,".txt"),quote = F,row.names = F,col.names = F )
  
}

extract.code <- rep("c",1000)
eth <- c("EUR","AFR","AMR","EAS")
temp <- 1
for(i in 1:length(eth)){
  for(j in 1:22){
    extract.code[temp] <- paste0("/data/zhangh24/software/qctool/qctool -g /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".combined.tag.gen -og /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_cau.combined.tag.gen -incl-positions /data/zhangh24/multi_ethnic/result/LD_simulation/causal_snp_pos_chr_",j,".txt")  
    temp <- temp+1
  }
}
extract.code <- extract.code[1:(temp-1)]
write(extract.code,file = "./code/LD_simulation/extract_causal_SNP.sh")

#there is bug for EUR chromsome 2 file by qctool
#manually extract SNPs using R
eth <- c("EUR","AFR","AMR","EAS")
i <- 1
j <- 2
library(data.table)
gen.infor <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/chr2.tag.info.txt")))
#the causal SNPs position
jdx <- which(cau.snp.infor$CHR==j)
snp.pos <-cau.snp.infor[jdx,]$position


idx <- which(gen.infor$V3%in%snp.pos)

gen.file <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".combined.tag.gen")
jdx <- which(cau.snp.infor$CHR==j)
snp.pos <-cau.snp.infor[jdx,]$position

result.list <- NULL
con <- file(gen.file)
total <- 0
open(con)
for(k in 1:nrow(gen.infor)){
  print(k)
  temp.line <- readLines(con,n=1)
  
  if(k%in%idx){
    result.list[[total+1]] <- temp.line
    total <- total+1
  }
}
close(con)

out.file <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_cau.combined.tag.gen")
sink(out.file)

for(k in 1:total){
  cat(result.list[[k]])
  cat("\n")
  
}

#cat(result.list[[total]])
sink()

#add the first NA column for the file

out.result <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_cau.combined.tag.gen"),header=F))
out.result <- cbind(NA,out.result)
dim(out.result)
write.table(out.result,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_cau.combined.tag.gen"),row.names = F,col.names = F,quote=F)






#merge all causal files into one causal SNPs file
merge.code <- rep("c",4)
eth <- c("EUR","AFR","AMR","EAS")
temp <- 1
for(i in 1:length(eth)){
  temp.code <- "cat"
  for(j in 1:22){
    temp.code <- paste0(temp.code," /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_cau.combined.tag.gen ")  
  }
  temp.code <- paste0(temp.code,
                      " > /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/cau.combined.tag.gen")
  merge.code[i] <- temp.code 
}
write(merge.code,file = 
        "./code/LD_simulation/merge_extract_SNP.sh")

