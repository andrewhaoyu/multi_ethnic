#get the duplicated SNP id for each population
eth <- c("EUR","AFR","AMR","EAS")
duplicated.id <- NULL
for(i in 1:length(eth)){
  for(j in 1:22){
    print(j)
    infor <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag.info.txt"),header=F))
    idx <- which(duplicated(infor[,2]))
    #EUR chr2 gen file get corrupted after line 697555
    #drop all of them from 697555 to 700000
    if(i==1&j==2){
      idx <- infor[697555:nrow(infor),2]
    }
    if(length(idx)!=0){
      duplicated.id <- c(duplicated.id,infor[idx,2]) 
    }
  }
  write.table(duplicated.id,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id"),row.names = F,col.names = F,quote=F)  
}
