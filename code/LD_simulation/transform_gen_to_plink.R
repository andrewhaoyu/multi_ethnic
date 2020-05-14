#transform gen file to plink file
#this way can make the run gwas and calculate prs faster since plink doens't need to reformat the gen to plink format again and again
code <- rep("c",1000)
temp <- 1
eth <- c("EUR","AFR","AMR","EAS")
for(i in 1:4){
  for(j in 1:22){
    code[temp] <- paste0("/data/zhangh24/software/plink2 --gen /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".plink.tag.gen --sample /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/sample.txt --make-bed --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag")
    
    
    temp <- temp+1
  }
}
code <-code[1:(temp-1)]
write.table(code,file="/data/zhangh24/multi_ethnic/code/LD_simulation/transform_gen_to_plink.sh",row.names = F,col.names = F,quote=F)



#reaload the fam file 
n.train = c(100000,15000,15000,15000)
for(i in 1:4){
  for(j in 1:22){
    fam <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag.fam"),header=F))
    load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))
    fam$V6[1:n.train[i]] <- y[1:n.train[i]]
    write.table(fam,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag.fam"),
                row.names = F,col.names = F,quote = F)
  }
}
