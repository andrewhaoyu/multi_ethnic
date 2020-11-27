#LD_clumping for different ethnic groups
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture
#r_ind represent r2 clumping threshold
#k_ind represent kbp clumping threshold
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
#r_ind = as.numeric(args[[6]])

library(data.table)
library(dplyr)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bed /lscratch/",sid,"/test/",eth[i],"_chr_all.bed"))
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bim /lscratch/",sid,"/test/",eth[i],"_chr_all.bim"))
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.fam /lscratch/",sid,"/test/",eth[i],"_chr_all.fam"))

system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bed ",temp.dir,eth[i],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bim ",temp.dir,eth[i],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.fam ",temp.dir,eth[i],"clump_ref_all_chr.fam"))

#update the summary results to make it work for plink clumping command

# for(i in 1:4){
#   summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))
#   colnames(summary)[1] <- "CHR"
#   assoc = summary %>%
#    # rename(SNP=ID,STAT=T_STAT,BP=POS) %>%
#     select(CHR,SNP,BP,A1,TEST,BETA,STAT,P)
#   write.table(assoc,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/assoc.out"),col.names = T,row.names = F,quote=F)
# 
# }
#only keep the SNPs with MAF > 0.05 in each ethnic group
setwd("/data/zhangh24/multi_ethnic/")
#load 1kg snp list
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
#idx <- which(snp.infor$position==752566)
#load mega snp list
snp.infor <- snp.infor.match
library(data.table)
library(dplyr)
MAF = snp.infor %>% 
  select(c(id,all_of(eth),rs_id)) %>% 
  rename(SNP=id)
# colnames(MAF)[1] <- "SNP"
# MAF_var = list('EUR','AFR','AMR','EAS',"SAS")
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
var_name = colnames(sum.data)
sum.data = left_join(sum.data,MAF,by="SNP")
library(tidyr)
mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(mega.list) = "rs_id"
sum.data.com <- inner_join(mega.list,sum.data,by="rs_id")


sum.data.MAF = sum.data.com %>% filter(get(eth[i])>=0.01&get(eth[i])<=0.99) %>% 
  select(all_of(var_name))
#idx.order <- order(sum.data.MAF$P)
#sum.data.MAF.new <- sum.data.MAF[idx.order,]
write.table(sum.data.MAF,file = paste0(temp.dir,eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
            ,col.names = T,row.names = F,quote=F) 

# dim(summary)
# head(summary)
r2_vec = c(0.01,0.05,0.1,0.2,0.5)
wc_base_vec = c(50,100,200,500)
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 0.5
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[wc_ind]
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
    #code <- rep("c",5*3*3)
    #system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[i],"clump_ref_all_chr --clump ",temp.dir,eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind))
    if(res==2){
      stop()
    }
    system(paste0("mv ",temp.dir,"LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped ",out.dir,eth[i],"/"))
    #system(paste0("rm ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist,".log"))
    
  }
}


