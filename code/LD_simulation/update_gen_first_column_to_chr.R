#plink2 require the first column of the gen file to be chromosome
#currently it's snpid
#we should update it using awk lanuage

code <- rep("c",1000)
eth <- c("EUR","AFR","AMR","EAS")
temp = 1
for(i in 1:length(eth)){
  for(j in 1:22){
    code[temp] <- paste0("awk '{$1=",j," ; print ;}' /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".combined.tag.gen > /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".plink.tag.gen")
    temp <- temp+1
  }
}
code = code[1:(temp-1)]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/update_gen.sh",quote = F,row.names = F,col.names = F)
