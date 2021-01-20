#LD_clumping for different ethnic groups combined with EUR
#LD_clumping for EUR
#load the p-value results
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
library(data.table)
library(dplyr)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/LD/'),showWarnings = FALSE)
temp.dir.LD <- paste0('/lscratch/',sid,'/test/LD/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#clumping for EUR SNPs
summary.eur <- as.data.frame(fread(paste0(out.dir,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
#only keep snpid and p-value for future combination
summary.eur.select = summary.eur %>% 
  select(SNP,peur)
#for(i in 2:4){
#summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))
summary <- as.data.frame(fread(paste0(out.dir,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
#find the shared SNPs between target ethnic group and EUR
#get the min p-value for between the target ethnic group and EUR for shared snp
summary.com <- left_join(summary,summary.eur.select,by="SNP")
#select the SNPs from EUR p-value
idx <- which(summary.com$peur<summary.com$P)
summary.com.EUR <- summary.com[idx,]

#load 1kg snp list
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
summary.com.EUR = left_join(summary.com.EUR,snp.infor,by="SNP")
mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(mega.list) = "rs_id"

#subset SNPs to mega SNPs
summary.match.EUR = inner_join(summary.com.EUR,mega.list,by="rs_id")
assoc = summary.match.EUR %>%
  select(CHR,SNP,BP,A1,TEST,BETA,STAT,peur) %>% 
  rename(P=peur)



write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[1],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)

system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bed ",temp.dir,eth[1],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bim ",temp.dir,eth[1],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.fam ",temp.dir,eth[1],"clump_ref_all_chr.fam"))
setwd("/data/zhangh24/multi_ethnic/")
r2_vec = c(0.01,0.05,0.1,0.2,0.5)
wc_base_vec = c(50,100)
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 0.5
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
    #code <- rep("c",5*3*3)
    #system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[1],"clump_ref_all_chr --clump ",temp.dir,eth[1],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[1],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind))
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
summary.com.tar <- summary.com[idx,]
summary.com.tar = left_join(summary.com.tar,snp.infor,by="SNP")

#subset SNPs to mega SNPs
summary.match.tar = inner_join(summary.com.tar,mega.list,by="rs_id")


assoc = summary.match.tar %>%
  select(CHR,SNP,BP,A1,TEST,BETA,STAT,P) 
#%>% 
#rename(P=peur)
write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bed ",temp.dir,eth[i],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bim ",temp.dir,eth[i],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.fam ",temp.dir,eth[i],"clump_ref_all_chr.fam"))


for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 0.5
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[i],"clump_ref_all_chr --clump ",temp.dir,eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[i],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind))
    if(res==2){
      stop()
    }
    
     }
}

for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
LD.EUR <- as.data.frame(fread(paste0(temp.dir,eth[1],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
LD.tar <- as.data.frame(fread(paste0(temp.dir,eth[i],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))

LD <- rbind(LD.EUR,LD.tar)
write.table(LD,file = paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped"),row.names = F,quote=F)
  }
  }
#LD.test <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".clumped")))
# summary.com = summary.com %>% 
#   mutate(p_update=pmin(P,peur,na.rm = T))
# write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
# 
# #update the summary results to make it work for plink clumping command
# 
# sid<-Sys.getenv('SLURM_JOB_ID')
# dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bed /lscratch/",sid,"/test/",eth[i],"_chr_all.bed"))
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bim /lscratch/",sid,"/test/",eth[i],"_chr_all.bim"))
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.fam /lscratch/",sid,"/test/",eth[i],"_chr_all.fam"))
# setwd("/data/zhangh24/multi_ethnic/")
# #summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary.out"),header=T))
# 
# #}
# 
# 
# dim(summary)
# head(summary)
# pthr = 0.5
# r2thr = 0.1
# kbpthr = 500
# res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/",eth[i],"_chr_all --clump /lscratch/",sid,"/test/",eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", out.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
# if(res==2){
#   stop()
# }
#system(paste0("more ",out.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".log"))
#system(paste0("rm ",out.dir,eth[i],"LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".log"))
# }
# #no need for rerunning in EUR
# code <- code[2:4]
# write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/LD_clumping_two_dim.sh",row.names=F,col.names=F,quote=F)
# #the p-value takes the min(AFR,EUR)
# #if the SNP only exist
# 


# data <- fread(paste0("/data/zhangh24/KG.plink/",eth[i],"/chr_all.bim"))
# 
# temp = seq(-0.01788854,0.01788854,by=0.0001)
# p.vec = mean(pnorm(-abs(sqrt(15000)*temp),lower.tail = T)/pnorm(-abs(sqrt(100000)*temp),lower.tail = T))
# 
