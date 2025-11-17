old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
for(i in 1:5){
  data <- as.data.frame(fread(paste0(old.out.dir,eth[i],"/duplicated.id"),header=F))
  data <- unique(data)
  write.table(data,file = paste0(old.out.dir,eth[i],"/duplicated.id"),row.names = F,col.names = F,quote = F)
  
}
