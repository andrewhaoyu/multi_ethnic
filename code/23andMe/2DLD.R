#LD_clumping for ARIC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
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
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bed ",temp.dir,eth[1],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bim ",temp.dir,eth[1],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.fam ",temp.dir,eth[1],"all_chr.fam"))
method = "TDLD"
# for(i in 1:2){
#   for(l in 1:3){
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
#out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")
#load EUR gwas summary statistics
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
#select the SNPs from EUR p-value
idx <- which(summary.com$peur<summary.com$P)
summary.com.EUR <- summary.com[idx,] %>% 
  select(CHR,SNP,BP,A1,BETA,peur) %>% 
  rename(P = peur)

#sum.data.assoc = sum.data.assoc[,c("CHR","SNP","BP","A1","BETA","P")]
#idx <- which(sum.data.assoc$SNP=="rs4970836")
write.table(summary.com.EUR,file = paste0(temp.dir,eth[1],"_all_chr_assoc.txt"),col.names = T,row.names = F,quote=F)
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bed ",temp.dir,eth[1],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bim ",temp.dir,eth[1],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.fam ",temp.dir,eth[1],"all_chr.fam"))
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
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
  }
}
system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))


idx <- which((summary.com$P<summary.com$peur)|is.na(summary.com$peur))
sum.com.tar <- summary.com[idx,] %>% 
  select(CHR,SNP,BP,A1,BETA,P) 
write.table(sum.com.tar,file = paste0(temp.dir,eth[i],"_all_chr_assoc.txt"),col.names = T,row.names = F,quote=F)
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bed ",temp.dir,eth[i],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bim ",temp.dir,eth[i],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.fam ",temp.dir,eth[i],"all_chr.fam"))


for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
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
    
  }
}


system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))



out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")

for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    LD.EUR <- as.data.frame(fread(paste0(temp.dir,eth[1],"_TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    LD.tar <- as.data.frame(fread(paste0(temp.dir,eth[i],"_TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    
    LD.EUR <- LD.EUR[,3,drop=F]
    LD.tar <- LD.tar[,3,drop=F]
    LD <- rbind(LD.EUR,LD.tar)
    write.table(LD,file = paste0(out.dir,"TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped"),row.names = F,quote=F)
  }
}


#   }
# }
