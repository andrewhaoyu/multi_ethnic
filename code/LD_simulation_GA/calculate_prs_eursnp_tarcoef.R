#select SNPs based on best EUR prs, use coefficients of target population
#i represent ethnic group
#l represent causal proportion
#m represent sample size
args = commandArgs(trailingOnly = T)
i_rep = as.numeric(args[[1]])
i = as.numeric(args[[2]])
j = as.numeric(args[[3]])
i1 = as.numeric(args[[4]])
setwd("/data/zhangh24/multi_ethnic/")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)

sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0("/lscratch/",sid,"/",eth[i],"/"),showWarnings = F)
temp.dir <- paste0("/lscratch/",sid,"/",eth[i],"/")
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.* ",temp.dir))
system(paste0("ls ",temp.dir))
library(dplyr)
library(data.table)
for(l in 1:3){
  for(m in 1:4){
    #n.snp.mat <- matrix(0,length(pthres),4)
    #load EUR r2 result
    load(paste0(out.dir,"LD.clump.result.GA_",i1,".rdata"))
    #keep the EUR results with sample size at 100,000 (m = 4)
    r2.mat <- as.data.frame(LD.result.list[[2]]) %>% 
      filter(eth.vec=="EUR"&
               m_vec==4&
               l_vec==l)
    
    #use target population regresssion coefficients
    #get the number of SNPs in different population
    #load EUR SNP
    k = which.max(r2.mat$r2.vec)
    LD <- as.data.frame(fread(paste0(out.dir,eth[1],"/LD_clump_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1,".clumped")))
    clump.snp <- LD[,3,drop=F]  
    #select the EUR snp
    sum.data <- as.data.frame(fread(paste0(out.dir,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
    colnames(sum.data)[2] <- "SNP"
    prs.all <- left_join(clump.snp,sum.data,by="SNP") 
    
    prs.file <- prs.all %>% filter(P<=pthres[k]) %>% 
      select(SNP)
    sum.data <- as.data.frame(fread(paste0(out.dir,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
    #combine the prs with the target population coefficients
    #idx <- which(prs.file$SNP%in%sum.data$SNP==F)
    #length(idx)
    prs.file.com <- left_join(prs.file,sum.data,by="SNP")
    prs.file <- prs.file.com %>% select(SNP,A1,BETA)
    #prs.file[14,]
    prs.file <- prs.file[complete.cases(prs.file),]
    write.table(prs.file,file = paste0("/lscratch/",sid,"/",eth[i],"/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --score /lscratch/",sid,'/',eth[i],"/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1, " no-sum no-mean-imputation --bfile ",temp.dir,"/chr",j,".tag --exclude ",out.dir,eth[i],"/duplicated.id  --out ",out.dir,eth[i],"/prs/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1))
    if(res==2){
      stop()
    }
    system(paste0("rm ",cur.dir,eth[i],"/prs/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".log"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".nosex"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".nopred"))
    system(paste0("rm /lscratch/",sid,"/",eth[i],"/prs_eursnp_tarcoef_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    gc()
  }
}
#n.snp.mat <- matrix(0,length(pthres),4)
#load EUR r2 result
# bim.file <- fread(paste0(temp.dir,"/chr",j,".tag.bim"),header=F)
# colnames(bim.file)[2] <- "SNP"
# new.file <- inner_join(prs.file,bim.file,by="SNP")
# dim(new.file)
#idx <- which(prs.file$SNPas.character(bim.file[,2]))
#length(idx)
system(paste0('rm -r /lscratch/',sid,'/',eth[i],'/'))