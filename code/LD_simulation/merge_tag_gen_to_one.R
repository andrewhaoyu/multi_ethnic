#merge all the 100 replicates tag gene files into one
#i1 represent ethnic groups
#i2 represent chromosome
eth <- c("EUR","AFR","AMR","EAS")
#tag <- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i1],"_chr",i2,".tag"),header=F)
code <- rep("c",1000)
temp <- 1
for(i1 in 1:1){
  for(i2 in 1:22){
    sub.code <- paste0("paste -d \" \" /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",i2,".tag.info.txt")
    for(i3 in 1:500){
      sub.code <- paste0(sub.code," /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",i2,"_",i3,".controls.tag.gen")
    }
    sub.code <- paste0(sub.code,
                       " > /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",i2,".combined.tag.gen")
    code[temp] <- sub.code
    temp <- temp+1    
  }

}
code <- code[1:(temp-1)]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/merge_tag_gen_to_one.sh",quote=F,
            row.names = F,col.names = F)



#check whether the files exist
i1 <- 1
filedir <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1])
files <- dir(filedir,pattern = ".tag.gen")

resubmit.file <- NULL

for(i2 in 1:22){
  for(i3 in 1:500){
 temp.file <-    paste0("chr",i2,"_",i3,".controls.tag.gen")
    if(temp.file%in%files==F){
      resubmit.file <- c(resubmit.file,temp.file)
    }
  }
}
  








