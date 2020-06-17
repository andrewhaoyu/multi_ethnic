#use plink2 to calculate prs
#load LD_clump_file
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(1E-10,1E-09,5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03)
#pthres <- c(1E-10,1E-09,5E-08,1E-07,2.5E-07,5E-07,7.5E-07,1E-06,2.5E-06,5E-06,7.5E-06,1E-05,2.5e-05,5E-05,7.5e-05,1E-04,2.5E-04,5E-04,7.5E-04,1E-03)
n.snp.mat <- cbind(expand.grid(pthres,pthres),0,0,0)
colnames(n.snp.mat) = c("pthres1","pthres2","nsnp")
#read summary level statistics from EUR
summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary.out"),header=T))
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>% 
  select(SNP,beta_eur,peur)
for(i in 2:length(eth)){
  #read LD clumped SNPs from EUR
  LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump_eurp.clumped")))
  clump.snp <- LD[,3,drop=F] 
  #read the target ethnic group summary level statistics
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  colnames(sum.data)[2] <- "SNP"
  #combine the target level summary stat with EUR
  summary.com <- left_join(sum.data,summary.eur.select,by="SNP")
  #combine the statistics with SNPs after clumping
  prs.all <- left_join(clump.snp,summary.com,by="SNP") 
  #save the prs information for easier filter
  save(prs.all,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump_two_dim.clumped_all_infor_eurp_tarref"))
  
}

for(i in 2:length(eth)){
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump_two_dim.clumped_all_infor_eurp_tarref"))
  #k for the p-value threshold on the target population
  #l for the p-value threshold on the eur population
  for(k in 1:length(pthres)){
    for(l in 1:length(pthres)){
      prs.file <- prs.all %>% filter(P<=pthres[k]|
                                       peur<=pthres[l]) %>% 
        select(SNP,A1,BETA)
      dim(prs.file)
      
      #       prs.file.new <- prs.all %>% filter(peur<=pthres[l]) %>% 
      #         select(SNP,A1,BETA)
      # dim(prs.file.new)      
      
      
      jdx = which(n.snp.mat[,1]==pthres[k]&n.snp.mat[,2]==pthres[l])
      n.snp.mat[jdx,i+1] = nrow(prs.file)
      write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_two_dim_eurp_tarref",k,"_",l),col.names = T,row.names = F,quote=F)  
    }
    
  }
  
}
write.csv(n.snp.mat,file =paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/n_snp_mat.csv"))

#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 2:length(eth)){
  for(j in 1:22){
    for(k in 1:length(pthres)){
      for(l in 1:length(pthres)){
        temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_two_dim_eurp_tarref",k,"_",l," no-sum no-mean-imputation  --allow-no-sex --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr_eurp_tarref",j,"_prs_",k,"_",l)
        code[temp] <- temp.code
        temp <- temp+1  
      }
    }
  }
  
}
code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs_two_dim_eurp_tarref.sh"),col.names = F,row.names = F,quote=F)







