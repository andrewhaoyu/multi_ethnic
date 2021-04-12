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
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bed ",temp.dir,eth[i],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bim ",temp.dir,eth[i],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.fam ",temp.dir,eth[i],"all_chr.fam"))

# for(i in 1:2){
#   for(l in 1:3){
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/PT/",eth[i],"/",trait[l],"/")
#load gwas summary statistics
sum.data = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
# write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
#             ,col.names = T,row.names = F,quote=F) 
#prepare association file for plink
sum.data.assoc = sum.data %>% 
  rename(SNP = rsid) %>% 
  select(CHR,SNP,BP,A1,BETA,P) 

#sum.data.assoc = sum.data.assoc[,c("CHR","SNP","BP","A1","BETA","P")]
#idx <- which(sum.data.assoc$SNP=="rs4970836")
write.table(sum.data.assoc,file = paste0(temp.dir,"all_chr_assoc.txt"),col.names = T,row.names = F,quote=F)

# dim(summary)
# head(summary)
pthr = 1
r2thr = 0.1
kbpthr = 500
eth <- c("EUR","AFR","AMR","EAS","SAS")
#cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#code <- rep("c",5*3*3)
#system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
res = system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,eth[i],"all_chr --clump ",temp.dir,"all_chr_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_clump"))
system(paste0("mv ",temp.dir,"LD_clump.clumped ",out.dir))
if(res==2){
  stop()
}

#   }
# }