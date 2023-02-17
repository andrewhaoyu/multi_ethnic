args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
setwd("/data/zhangh24/multi_ethnic/")
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
method = "TDLD"
#for(i in 1:length(eth)){
 # for(l in 1:length(trait)){
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")
    #load gwas summary statistics
    sum.eur = as.data.frame(fread(paste0(data.dir,eth[1],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
    # write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
    #             ,col.names = T,row.names = F,quote=F) 
    #prepare association file for plink
    sum.eur.select = sum.eur %>% 
      rename(SNP = rsid,
             peur = P) %>% 
      select(SNP,peur) 
    
    
    #load target ethnic group data
    sum.tar = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
    
    sum.data.assoc = sum.tar %>% 
      rename(SNP = rsid) %>% 
      select(CHR,SNP,BP,A1,BETA,P) 
    
    #get the min p-value for between the target ethnic group and EUR for shared snp
    summary.com <- left_join(sum.data.assoc,sum.eur.select,by="SNP")
    temp = 1
    out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
    for(r_ind in 1:length(r2_vec)){
      wc_vec = wc_base_vec/r2_vec[r_ind]
      for(w_ind in 1:length(wc_vec)){
    LD <- as.data.frame(fread(paste0(out.dir,"TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD
    prs.all <- left_join(clump.snp,summary.com)
    
    
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        prs.temp <- prs.all %>% 
          filter((P<=pthres[k1]|
                    peur<=pthres[k2])) %>% 
          select(SNP,BETA)
        if(nrow(prs.temp)>0){
          write.table(prs.temp,file = paste0(out.dir.prs,temp,"_",method,"_rind_",r_ind,"_wcind_",w_ind,"_ptar_",k1,"_peur_",k2),row.names = F,col.names = T,quote=F)  
          temp = temp + 1
        }
        
      }
      
    }
      }
    }
 # }
#}
