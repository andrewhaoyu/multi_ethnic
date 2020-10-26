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
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#update the summary results to make it work for plink clumping command
eth <- c("EUR","AFR","AMR","EAS","SAS")
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
eth <- c("EUR","AFR","AMR","EAS","SAS")
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bed /lscratch/",sid,"/test/",eth[i],"_chr_all.bed"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bim /lscratch/",sid,"/test/",eth[i],"_chr_all.bim"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.fam /lscratch/",sid,"/test/",eth[i],"_chr_all.fam"))
setwd("/data/zhangh24/multi_ethnic/")
#summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary.out"),header=T))
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
summary.com = summary.com %>% 
  mutate(p_update=pmin(P,peur,na.rm = T))
assoc = summary.com %>%
  # rename(SNP=ID,STAT=T_STAT,BP=POS) %>%
  select(CHR,SNP,BP,A1,TEST,BETA,STAT,p_update) %>% 
  rename(P=p_update)
write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)

#}


dim(summary)
head(summary)
pthr = 0.5
r2thr = 0.1
kbpthr = 500
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/",eth[i],"_chr_all --clump /lscratch/",sid,"/test/",eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", out.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
if(res==2){
  stop()
}
#system(paste0("more ",out.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".log"))
system(paste0("rm ",out.dir,eth[i],"LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".log"))
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
