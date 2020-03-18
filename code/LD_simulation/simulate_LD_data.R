#simulate genotype data by chromosome and by ethnic groups
#EUR have 120*k genotype data
#AFR, AMR, EAS, SAS have 18k genotype data

#since hapgen2 required at least one disease SNP, we need to know the position information for one SNP

#separate EUR and other three ethnic groups
#since EUR requires more jobs

#simulate data for AFR, AMR, EAS
library(data.table)
eth <- c("AFR","AMR","EAS")

code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(180,180,180)
tag_all <- list()
for(i in 1:length(eth)){
  tag_temp <- list()
  for(j in 1:22){
    tag_temp[[j]] <- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i],"_chr",j,".tag"),header=F)
  }
  tag_all[[i]] <- tag_temp
  
}

#for put chr and ethnic groups in inner loop to aviod the same start time for hapgen2
#hapgen2 use the start time to set random seed
for(k in 1:100){
  for(j in 1:22){
    
    for(i in 1:length(eth)){    
      tag <- tag_all[[i]][[j]]
      code[temp] <- paste0("/data/zhangh24/software/hapgen2 -m /data/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt -l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend -h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap -o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -no_haps_output"
      )  
      temp = temp+1   
    }
    
    
  }
  
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data_other_eth.sh"),row.names = F,col.names = F,quote=F)


#simulate data for EUR

eth <- c("EUR")

code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(120)
tag_all <- list()
for(i in 1:length(eth)){
  tag_temp <- list()
  for(j in 1:22){
    tag_temp[[j]] <- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i],"_chr",j,".tag"),header=F)
  }
  tag_all[[i]] <- tag_temp
  
}

#for put chr and ethnic groups in inner loop to aviod the same start time for hapgen2
#hapgen2 use the start time to set random seed
for(k in 1:1000){
  for(j in 1:22){
    
    for(i in 1:length(eth)){    
      tag <- tag_all[[i]][[j]]
      code[temp] <- paste0("/data/zhangh24/software/hapgen2 -m /data/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt -l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend -h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap -o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -no_haps_output"
      )  
      temp = temp+1   
    }
    
    
  }
  
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_EUR.sh"),row.names = F,col.names = F,quote=F)








code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(120,180,180,180,180)

for(j in 22){
  tag <- read.table(paste0("/spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag"),header=F)
  for(k in 1:100){
    for(i in 1){
      code[temp] <- paste0("TZ='EDT+",format(Sys.time()+temp*10,"%H:%M:%S"),"' /data/zhangh24/software/hapgen2 -m /spin1/users/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt -l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend -h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap -o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -t /spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag -no_haps_output"
      )  
      temp = temp+1   
    }
    
    
  }
  
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data.sh"),row.names = F,col.names = F,quote=F)

