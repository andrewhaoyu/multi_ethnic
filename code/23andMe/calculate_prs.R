#args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

# i = as.numeric(args[[1]])
# l = as.numeric(args[[2]])
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
method = "PT"
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"

for(i in 1:length(eth)){
  for(l in 1:length(trait)){
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/PT/",eth[i],"/",trait[l],"/")
    #load gwas summary statistics
    sum.data = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
    # write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
    #             ,col.names = T,row.names = F,quote=F) 
    #prepare association file for plink
    sum.data.assoc = sum.data %>% 
      rename(SNP = rsid) %>% 
      select(CHR,SNP,BP,A1,BETA,P) 
    LD <- as.data.frame(fread(paste0(out.dir,"LD_clump.clumped")))
    clump.snp <- LD[,3,drop=F]
    prs.all <- left_join(clump.snp,sum.data.assoc)
    temp = 1
    out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PT/",eth[i],"/",trait[l],"/")
    
    for(k in 1:length(pthres)){
      prs.temp <- prs.all %>% 
        filter(P<=pthres[k]) %>% 
        select(SNP,BETA)
      if(nrow(prs.temp)>0){
        write.table(prs.temp,file = paste0(out.dir.prs,temp,"_",method,"_pvalue_",k),row.names = F,col.names = T,quote=F)  
        temp = temp + 1
      }
      
    }
    
  }
}
