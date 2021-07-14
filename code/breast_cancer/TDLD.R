args = commandArgs(trailingOnly = T)
r_ind = as.numeric(args[[1]])
w_ind = as.numeric(args[[2]])
l = 1
i = 2
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
method = "TDLD"
#load EUR gwas summary statistics

setwd("/data/zhangh24/multi_ethnic/data/")
load("./AABC_data/BC_EUR_overall_mega_aligned.rdata")
sum.eur = sum.data
trait = c("overall","erpos","erneg")
#load target ethnic group data
sum.eur.select = sum.eur %>% 
  rename(SNP=V1,
         peur = p_eur) %>% 
  select(SNP,peur) 
load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS_mega.rdata"))
sum.tar = sum.data
sum.data.assoc = sum.tar %>% 
  rename(SNP = V1,
         A1 = Effect_allele,
         BETA = Effect,
         BP = POS) %>% 
  select(CHR,SNP,BP,A1,BETA,P,ID) 

summary.com <- left_join(sum.data.assoc,sum.eur.select,by="SNP")
#select the SNPs from EUR p-value
idx <- which(summary.com$peur<summary.com$P)
summary.com.EUR <- summary.com[idx,] %>% 
  select(CHR,SNP,BP,A1,BETA,peur) %>% 
  rename(P = peur)
write.table(summary.com.EUR,file = paste0(temp.dir,eth[1],"_all_chr_assoc.txt"),col.names = T,row.names = F,quote=F)
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bed ",temp.dir,eth[1],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bim ",temp.dir,eth[1],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.fam ",temp.dir,eth[1],"all_chr.fam"))
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
# for(r_ind in 1:length(r2_vec)){
 wc_vec = wc_base_vec/r2_vec[r_ind]
#   for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 1
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    
    #code <- rep("c",5*3*3)
    #system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[1],"all_chr --clump ",temp.dir,eth[1],"_all_chr_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[1],"_TDLD_rind_",r_ind,"_wcind_",w_ind))
    if(res==2){
      stop()
    }
    #system(paste0("mv ",temp.dir,"LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped ",out.dir,eth[i],"/"))
#   }
# }
system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))

idx <- which((summary.com$P<summary.com$peur)|is.na(summary.com$peur))
#Ghana data used ID instead of rsid
sum.com.tar <- summary.com[idx,] %>% 
  select(CHR,ID,BP,A1,BETA,P) %>% 
  rename(SNP =ID)
write.table(sum.com.tar,file = paste0(temp.dir,eth[i],"_all_chr_assoc.txt"),col.names = T,row.names = F,quote=F)
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bed ",temp.dir,eth[i],"all_chr.bed"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bim ",temp.dir,eth[i],"all_chr.bim"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.fam ",temp.dir,eth[i],"all_chr.fam"))

#for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
 # for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 1
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    #eth <- c("EUR","AFR","AMR","EAS","SAS")
    #out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[i],"all_chr --clump ",temp.dir,eth[i],"_all_chr_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[i],"_TDLD_rind_",r_ind,"_wcind_",w_ind))
    if(res==2){
      stop()
    }
    
#   }
# }

system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))



out.dir = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/clump_result/")



#for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
 # for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    LD.EUR <- as.data.frame(fread(paste0(temp.dir,eth[1],"_TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    LD.tar <- as.data.frame(fread(paste0(temp.dir,eth[i],"_TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    
    LD.EUR <- LD.EUR[,3,drop=F]
    LD.tar <- LD.tar[,3,drop=F]
    LD <- rbind(LD.EUR,LD.tar)
    write.table(LD,file = paste0(out.dir,"TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped"),row.names = F,quote=F)
  #}
#}
