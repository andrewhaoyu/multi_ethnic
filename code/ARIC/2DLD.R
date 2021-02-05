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
    #clumping based on EUR
    sum.data.assoc.EUR <- summary.com %>% 
      filter(peur<PVAL & CHR ==j) %>% 
      mutate(BP=POS, SNP = SNP_ID,
             A1 = REF,
        P = peur) %>% 
      select(CHR,SNP,BP,A1,BETA,P)
    
    #clump based on AFR
    sum.data.assoc.tar <- summary.com %>% 
      filter(((PVAL<peur)|is.na(peur)) & CHR ==j) %>% 
      mutate(BP=POS, SNP = SNP_ID,
             A1 = REF,
             P = PVAL) %>% 
      select(CHR,SNP,BP,A1,BETA,P)
    
    #idx <- which(sum.data.assoc$SNP=="rs4970836")
    write.table(sum.data.assoc.tar,file = paste0(temp.dir,"2DLD_chr_",j,"_assoc_tar.txt"),col.names = T,row.names = F,quote=F)
    
      
    #idx <- which(sum.data.assoc$SNP=="rs4970836")
    write.table(sum.data.assoc.EUR,file = paste0(temp.dir,"2DLD_chr_",j,"_assoc_EUR.txt"),col.names = T,row.names = F,quote=F)
    r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
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
    res = system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --bfile ",data.dir,trait[1],"/",eth[1],"/geno/mega/ref_chr",j," --clump ",temp.dir,"2DLD_chr_",j,"_assoc_EUR.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,"_EUR"))
    res = system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/ref_chr",j," --clump ",temp.dir,"2DLD_chr_",j,"_assoc_tar.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,"_tar"))
    #system(paste0("mv ",temp.dir,"2DLD_clump_chr_",j,".clumped ",out.dir))
    if(res==2){
      stop()
    }
    
    clump.eur.snp <- as.data.frame(fread(paste0(temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,"_EUR.clumped"))) %>% 
      select(SNP)
    
    clump.tar.SNP <- clump.eur <- as.data.frame(fread(paste0(temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,"_tar.clumped"))) %>% 
      select(SNP)
    clump.snp <- rbind(clump.eur.snp,clump.tar.SNP)
    write.table(clump.snp, file = paste0(temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,".clumped"))
    system(paste0("rm ",temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,"_EUR.clumped"))
    system(paste0("rm ",temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,"_tar.clumped"))
      }
    }
    
  }
  



    
  
