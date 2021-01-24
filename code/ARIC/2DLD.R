#LD_clumping for ARIC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

#i = as.numeric(args[[1]])
#l = as.numeric(args[[2]])
j = as.numeric(args[[1]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
i = 2
#for(i in 1:2){
  for(l in 1:3){
    setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    #load gwas data for EUR SNPs
    sum.eur = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[1],"/sumdata/training-GWAS-formatted.txt")))
    sum.eur = sum.eur %>% rename(peur = PVAL) %>% 
      select(SNP_ID,peur)
    
    #load target gwas summary statistics
    sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    #find the shared SNPs between target ethnic group and EUR
    #get the min p-value for between the target ethnic group and EUR for shared snp
    summary.com <- left_join(sum.data,sum.eur,by="SNP_ID")
    
    P = summary.com$PVAL
    idx <- which(summary.com$peur<summary.com$PVAL)    
    P[idx] = summary.com$peur[idx]
    summary.com$P = P
    
    sum.data.assoc = summary.com %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF) %>% 
      filter(CHR==j) 
    
    sum.data.assoc = sum.data.assoc[,c("CHR","SNP","BP","A1","BETA","P")]
    #idx <- which(sum.data.assoc$SNP=="rs4970836")
    write.table(sum.data.assoc,file = paste0(temp.dir,"2DLD_chr_",j,"_assoc.txt"),col.names = T,row.names = F,quote=F)
    r2_vec = c(0.01,0.05,0.1,0.2,0.5)
    wc_base_vec = c(50,100)
    for(r_ind in 1:length(r2_vec)){
      wc_vec = wc_base_vec/r2_vec[r_ind]
      for(w_ind in 1:length(wc_vec)){
        print(c(r_ind,w_ind))
        pthr = 0.5
        r2thr = r2_vec[r_ind]
        kbpthr = wc_vec[w_ind]
        
    #cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
    #code <- rep("c",5*3*3)
    #system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
    res = system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/ref_chr",j," --clump ",temp.dir,"2DLD_chr_",j,"_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind))
    #system(paste0("mv ",temp.dir,"2DLD_clump_chr_",j,".clumped ",out.dir))
    if(res==2){
      stop()
    }
    
      }
    }
  }
#}
#system(paste0("rm ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".log"))

# data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/geno/mega/ref_chr",j,".bim")))
# idx <- which(data$V2=="rs78444298")


# eth <- c("EUR","AFR","AMR","EAS","SAS")
# trait = c("eGFRcr","ACR","urate")
# out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/")
# for(l in 1:3){
#   system(paste0("mkdir ",out.dir,trait[l]))
#   for(i in 1:2){
#     system(paste0("mkdir ",out.dir,trait[l],"/",eth[i],"/"))
#   }
#}
#quality check


# all.bim.list = list()
# 
# for(j in 1:22){
#   file = paste0(data.dir,trait[l],"/",eth[i],"/geno/mega/ref_chr",j,".bim")
#   temp.bim <- as.data.frame(fread(file))
#   all.bim.list[[j]] = temp.bim
# }
# 
# all.bim = rbindlist(all.bim.list)
# 
# 
# sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
# sum.data.all = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-all.txt")))
# 
# head(sum.data.all)
# 
# sum.data.all = sum.data.all %>% 
#   mutate(SNP_ID = MARKER) %>% 
#   select(SNP_ID,EFFECTALLELE,OTHERALLELE,BETA,EAF)
# 
# sum.data.com = left_join(sum.data,sum.data.all,by="SNP_ID")
# head(sum.data.com)
# 
# sum.data.com.sig = sum.data.com %>% 
#   filter(PVAL<5E-08)
# 
# idx <- which(sum.data.com.sig$SNP_ID%in%all.bim$V2==F)
# save(sum.data.com.sig[idx,],file=
# 
#        
#        range(ifelse(sum.data.com.sig[idx,]$EAF>0.5,1-sum.data.com.sig[idx,]$EAF,sum.data.com.sig[idx,]$EAF))
#        
# hist(sum.data.com.sig[idx,]$EAF)
# #       temp = temp +1 
# #     }
# #   }
# #   
# #   
# # }
# # write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation_large/LD_clumping.sh",row.names=F,col.names=F,quote=F)
# #the p-value takes the min(AFR,EUR)
# #if the SNP only exist
# 
# 
# 
# # 
# # sum.data = sum.data %>% 
# #   mutate(chr.pos = paste0(CHR,":",POS))
# #load snp information file
#  load("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/snp.infor.rdata")
#  idx <- which(snp.infor$position==)
# # snp.infor  = snp.infor %>% 
# #   mutate(chr.pos = paste0(CHR,":",position)) %>% 
# #   select(id,a0,a1,AFR,EUR,chr.pos)
# # 
# # sum.data.com = left_join(sum.data,snp.infor,by="chr.pos")
# # head(sum.data.com)
# # range(sum.data.com$EUR,na.rm=T)
# # range(sum.data.com$AFR,na.rm= T)
# # idx <- which(sum.data.com$EUR==0)
# # sum.data.com[idx,]
# # #MAF = snp.infor[,..eth]
# # MAF = cbind(snp.infor$id,MAF)
# # colnames(MAF)[1] <- "SNP"
# # MAF_var = list('EUR','AFR','AMR','EAS',"SAS")
# # sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
# # var_name = colnames(sum.data)
# # sum.data.com = left_join(sum.data,MAF,bY="SNP")
# # idx <- which(sum.data.com$EAS>=0.05&sum.data.com$EAS<=0.95)
# # idx2 <- which(sum.data.com$P<=5e-08)
# # sum.data.com.temp <- sum.data.com[idx2,]
# # sum.data.MAF = sum.data.com %>% filter(get(eth[i])>=0.01&get(eth[i])<=0.99) %>% 
# #   select(all_of(var_name))
# #idx.order <- order(sum.data.MAF$P)
# #sum.data.MAF.new <- sum.data.MAF[idx.order,]
# 
# 
# idx <- which(sum.data.com$SNP_ID=="rs78444298")
# sum.data.com[idx,]
