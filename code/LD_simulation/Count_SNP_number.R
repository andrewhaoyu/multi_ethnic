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
