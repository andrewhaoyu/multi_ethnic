#LD_clumping for different ethnic groups
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
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
library(data.table)
library(dplyr)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
eth <- c("EUR","AFR","AMR","EAS","SAS")
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bed /lscratch/",sid,"/test/",eth[i],"_chr_all.bed"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bim /lscratch/",sid,"/test/",eth[i],"_chr_all.bim"))
system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.fam /lscratch/",sid,"/test/",eth[i],"_chr_all.fam"))

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
load("./result/LD_simulation_new/snp.infor.rdata")
MAF = snp.infor[,..eth]
MAF = cbind(snp.infor$id,MAF)
colnames(MAF)[1] <- "SNP"
MAF_var = list('EUR','AFR','AMR','EAS',"SAS")
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
var_name = colnames(sum.data)
sum.data.com = left_join(sum.data,MAF,bY="SNP")
# idx <- which(sum.data.com$EAS>=0.05&sum.data.com$EAS<=0.95)
# idx2 <- which(sum.data.com$P<=5e-08)
# sum.data.com.temp <- sum.data.com[idx2,]
sum.data.MAF = sum.data.com %>% filter(get(eth[i])>=0.01&get(eth[i])<=0.99) %>% 
  select(all_of(var_name))
#idx.order <- order(sum.data.MAF$P)
#sum.data.MAF.new <- sum.data.MAF[idx.order,]
write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
            ,col.names = T,row.names = F,quote=F) 

# dim(summary)
# head(summary)
pthr = 0.5
r2thr = 0.1
kbpthr = 500
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#code <- rep("c",5*3*3)
#system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/",eth[i],"_chr_all --clump /lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
system(paste0("rm ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".log"))
#       temp = temp +1 
#     }
#   }
#   
#   
# }
# write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation_large/LD_clumping.sh",row.names=F,col.names=F,quote=F)
#the p-value takes the min(AFR,EUR)
#if the SNP only exist


