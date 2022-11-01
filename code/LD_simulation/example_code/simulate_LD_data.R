#simulate genotype data by chromosome and by ethnic groups

#since hapgen2 required at least one disease SNP, we need to know the position information for one SNP

#simulate data for EUR,AFR, AMR, EAS, SAS
#-tag command allows you to subset SNPs to a subset
#but I notice that the -tag command makes hapgen2 really slow
#so I just simulate all the SNPs and then subset the tag snps by myself
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")

code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(180,180,180,180,180)
#tag SNPs are the SNPs we need
tag_all <- list()
for(i in 1:length(eth)){
  tag_temp <- list()
  for(j in 1:22){
    tag_temp[[j]] <- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i],"_chr",j,".tag"),header=F)
  }
  tag_all[[i]] <- tag_temp
  
}

#for put chr and ethnic groups in inner loop to avoid the same start time for hapgen2
#hapgen2 use the start time to set random seed
for(k in 1:100){ #k is the number of replicate
  for(j in 1:22){ #j is the chromosome number
     
    for(i in 1:5){    #i is the eth (EUR, AFR, AMR, EAS, SAS)
      tag <- tag_all[[i]][[j]]
      code[temp] <- paste0("/data/zhangh24/software/hapgen2 ",
                           "-m /data/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt ",
                           "-l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend ",
                           "-h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap ",
                           "-o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," ",
                           "-n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -no_haps_output"
      )  
      temp = temp+1   
    }
    
    
  }
  
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data_other_eth.sh"),row.names = F,col.names = F,quote=F)






#alternative way is to add the tag flag to hapgen2, but it's very slow
code <- rep("c",length(eth)*10000)
temp <- 1
n <- c(180,180,180,180,180)
for(k in 1:100){
for(j in 22){
  for(i in 1){
  tag <- read.table(paste0("/spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag"),header=F)
  
    
      code[temp] <- paste0("/data/zhangh24/software/hapgen2 ",
      "-m /spin1/users/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt ",
      "-l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend ",
      "-h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap ",
      "-o /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 ",
      "-t /spin1/users/zhangh24/KG.impute2/tag/chr",j,".tag -no_haps_output")
      temp = temp+1   
    }
    
    
  }
  
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data.sh"),row.names = F,col.names = F,quote=F)

