#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
k = as.numeric(args[[4]])
j = as.nunmeric(args[[5]])
#j = as.numeric(args[[3]])
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
#n.snp.mat <- matrix(0,length(pthres),4)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
setwd("/data/zhangh24/multi_ethnic/")

    LD <- as.data.frame(fread(paste0(cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,".clumped")))
    clump.snp <- LD[,3,drop=F]  
    sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m)))  
    colnames(sum.data)[2] <- "SNP"
    prs.all <- left_join(clump.snp,sum.data,by="SNP") 
    #for(k in 1:length(pthres)){
    prs.file <- prs.all %>% filter(P<=pthres[k]) %>% 
      select(SNP,A1,BETA)
    write.table(prs.file,file = paste0(cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m),col.names = T,row.names = F,quote=F)
    #}
    #for(j in 1:22){
    if(nrow(prs.file)>0){
      system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
    }
    
    #}
    
 







