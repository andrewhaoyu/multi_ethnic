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
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bed /lscratch/",sid,"/test/",eth[i],"_chr_all.bed"))
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.bim /lscratch/",sid,"/test/",eth[i],"_chr_all.bim"))
# system(paste0("cp /data/zhangh24/KG.plink/",eth[i],"/chr_all.fam /lscratch/",sid,"/test/",eth[i],"_chr_all.fam"))

system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.bed /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag.bed"))
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.bim /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag.bim"))
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.fam /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag.fam"))

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
load("./result/LD_simulation_new/snp.infor.rdata")
#idx <- which(snp.infor$position==752566)
#load mega snp list

library(data.table)
library(dplyr)
MAF = snp.infor[,..eth]
MAF = cbind(snp.infor$id,MAF)
colnames(MAF)[1] <- "SNP"
MAF_var = list('EUR','AFR','AMR','EAS',"SAS")
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
var_name = colnames(sum.data)
library(tidyr)
sum.data = sum.data %>% 
  separate(col=SNP,into=c("rs_id","position_new","allele1","allele2"),sep=":",remove=F)
jdx <- which(sum.data$CHR==1&sum.data$BP==8876162)


mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(mega.list) = "rs_id"
sum.data.match <- inner_join(mega.list,sum.data,by="rs_id")

idx <- which(mega.list$V1%in%sum.data$rs_id==F)
head(idx)
mega.list[sample(idx,10),]

library(tidyr)
sum.data.com = left_join(sum.data.match,MAF,bY="SNP") 


# idx <- which(sum.data.com$EAS>=0.05&sum.data.com$EAS<=0.95)
# idx2 <- which(sum.data.com$P<=5e-08)
# sum.data.com.temp <- sum.data.com[idx2,]
sum.data.MAF = sum.data.com %>% filter(get(eth[i])>=0.01&get(eth[i])<=0.99) %>% 
  select(all_of(var_name))
#idx.order <- order(sum.data.MAF$P)
#sum.data.MAF.new <- sum.data.MAF[idx.order,]
write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist,".out")
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
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/",eth[i],"_chr_all --clump /lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist))
if(res==2){
  stop()
}
system(paste0("rm ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist,".log"))
#       temp = temp +1 
#     }
#   }
#   
#   
# }
# write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation_large/LD_clumping.sh",row.names=F,col.names=F,quote=F)
#the p-value takes the min(AFR,EUR)
#if the SNP only exist







data <- fread2("/gpfs/gsfs11/users/zhangh24/multi_ethnic/result/LD_simulation_new/EUR/temp/chr_13.txt",header=T)








