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

assoc = summary.com.EUR %>%
  select(CHR,SNP,BP,A1,TEST,BETA,STAT,peur) %>% 
  rename(P=peur)
write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[1],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)

system(paste0("cp /data/zhangh24/KG.plink/",eth[1],"/chr_all.bed /lscratch/",sid,"/test/",eth[1],"_chr_all.bed"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[1],"/chr_all.bim /lscratch/",sid,"/test/",eth[1],"_chr_all.bim"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[1],"/chr_all.fam /lscratch/",sid,"/test/",eth[1],"_chr_all.fam"))
setwd("/data/zhangh24/multi_ethnic/")
temp.dir <- paste0("/lscratch/",sid,"/test/")
pthr = 0.5
r2thr = 0.1
kbpthr = 500
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[1],"_chr_all --clump ",temp.dir,eth[1],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[1],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
if(res==2){
  stop()
}
LD.EUR <- as.data.frame(fread(paste0(temp.dir,eth[1],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".clumped")))

#system(paste0("ls ",temp.dir))
#system(paste0("more ",temp.dir,"AFR_LD_clump_two_dim_rho_1_size_1_rep_1_GA_1.log"))
system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))

#select the SNPs from target ethnic p-value
idx <- which((summary.com$P<summary.com$peur)|is.na(summary.com$peur))
summary.com.tar <- summary.com[idx,]

assoc = summary.com.tar %>%
  select(CHR,SNP,BP,A1,TEST,BETA,STAT,P) 
#%>% 
  #rename(P=peur)
write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bed /lscratch/",sid,"/test/",eth[i],"_chr_all.bed"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bim /lscratch/",sid,"/test/",eth[i],"_chr_all.bim"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.fam /lscratch/",sid,"/test/",eth[i],"_chr_all.fam"))
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,eth[i],"_chr_all --clump ",temp.dir,eth[i],"_assoc.out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[i],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
if(res==2){
  stop()
}
LD.tar <- as.data.frame(fread(paste0(temp.dir,eth[i],"_LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".clumped")))

LD <- rbind(LD.EUR,LD.tar)
write.table(LD,file = paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".clumped"),row.names = F,quote=F)

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
