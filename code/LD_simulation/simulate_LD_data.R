#simulate genotype data by chromosome and by ethnic groups
#EUR have 120*k genotype data
#AFR, AMR, EAS, SAS have 18k genotype data

#since hapgen2 required at least one disease SNP, we need to know the position information for one SNP

library(data.table)

eth <- c("EUR","AFR","AMR","EAS","SAS")

code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(120,180,180,180,180)
  for(j in 1:22){
    tag <- read.table(paste0("/spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag"),header=F)
  for(k in 1:100){
    for(i in 1:length(eth)){
      code[temp] <- paste0("/data/zhangh24/software/hapgen2 -m /spin1/users/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt -l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend -h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap -o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -t /spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag -no_haps_output"
      )  
      temp = temp+1   
       }
  
 
  }
  
  
  }


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data.sh"),row.names = F,col.names = F,quote=F)


code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(120,180,180,180,180)
for(j in 22){
  tag <- read.table(paste0("/spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag"),header=F)
  for(k in 1:100){
    for(i in 1){
      code[temp] <- paste0("/data/zhangh24/software/hapgen2 -m /spin1/users/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt -l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend -h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap -o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -t /spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag -no_haps_output"
      )  
      temp = temp+1   
    }
    
    
  }
  
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data.sh"),row.names = F,col.names = F,quote=F)

