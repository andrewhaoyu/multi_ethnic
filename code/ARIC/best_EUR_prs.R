args = commandArgs(trailingOnly = T)
j = as.numeric(args[[1]])
library(dplyr)
library(data.table)

pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
  
for(l in 1:3){
    #load best EUR SNPs
    
    i = 2
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    out.dir.eur = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[1],"/")
    #for(j in 1:22){
    LD <- as.data.frame(fread(paste0(out.dir.eur,"/LD_clump_chr_",j,".clumped")))
    clump.snp <- LD[,3,drop=F] 
    
    sum.eur = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[1],"/sumdata/training-GWAS-formatted.txt")))
    sum.data = sum.eur %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL)
    prs.clump <- left_join(clump.snp,sum.data,by="SNP")
    
    #load LD clumping results
    load(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rdata"))
    CT.result = ARIC.result.CT[[2]] %>% 
      filter(triat==all_of(trait[l])&
               eth==all_of(eth[i]))
    idx.pcut <- which.max(CT.result$r2.vec.test.prs)
    prs.file = prs.clump %>% filter(P<=pthres[idx.pcut]) %>% 
      select(SNP,A1,BETA)
    #best eur SNP with EUR coefficients
    #write.table(prs.file,file = paste0(temp.dir,"best_eur_prs_coeff_chr_",j),col.names = T,row.names = F,quote=F)
    prs.eur = prs.file  
    
    #best EUR SNP with target coefficients
    sum.tar = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    sum.tar = sum.tar %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL)
    
    best.snp.id = prs.file %>% select(SNP)
    prs.file = left_join(best.snp.id,sum.tar,by="SNP") %>% 
      select(SNP,A1,BETA)
    prs.tar = prs.file
    #write.table(prs.file,file = paste0(temp.dir,"best_eur_tarcoef_prs_coeff_chr_",j),col.names = T,row.names = F,quote=F)
    #res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --threads 2 --score ",temp.dir,"best_eur_tarcoef_prs_coeff_chr_",j," header no-sum no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"best_eur_tarcoef_prs_chr_",j))
    
    #best EUR SNP with EB coefficients
    #align EUR and tar
    #use all best EUR SNPs to estimate prior
    clump.snp <- as.data.frame(fread(paste0(out.dir.eur,"/LD_clump.clumped")))
    sum.eur = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[1],"/sumdata/training-GWAS-formatted.txt")))
    sum.eur = sum.eur %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL)
    prs.clump <- left_join(clump.snp,sum.eur,by="SNP")
    best.snp = prs.clump %>% filter(P<=pthres[idx.pcut]) %>% 
      rename(A1_EUR=A1,BETA_EUR= BETA,
                 SE_EUR = SE) %>% 
      select(SNP,A1_EUR,BETA_EUR,SE_EUR,CHR)
    
    
    
    sum.tar = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    sum.tar = sum.tar %>% 
      mutate(BP=POS,SNP = SNP_ID,
             P = PVAL,
             BETA_tar = BETA,
             SE_tar = SE,
             A1_tar = REF) %>% 
      select(SNP,A1_tar,BETA_tar,SE_tar)
    
    sum.match <- left_join(best.snp,sum.tar,by="SNP")
    
    idx <- which(sum.match$A1_tar!=
                   sum.match$A1_EUR)
    if(length(idx)>0){
      sum.match$BETA_EUR[idx]=-sum.match$BETA_EUR[idx]
    }
    sum.match = sum.match %>% 
      mutate(z_stat_eur = BETA_EUR/SE_EUR,
             z_stat_tar = BETA_tar/SE_tar)
    prior.sigma = cov(cbind(sum.match$z_stat_tar,
                            sum.match$z_stat_eur),use="complete.obs")-diag(2)
    post.sigma = solve(solve(prior.sigma)+diag(2))
    z_stat_tar = sum.match$z_stat_tar
    z_stat_eur = sum.match$z_stat_eur
    z_mat = cbind(z_stat_tar,z_stat_eur)
    
    z_post = z_mat%*%post.sigma
    post_beta_tar = sum.match$BETA_tar
    idx <- which(!is.na(sum.match$z_stat_tar))
    post_beta_tar[idx] = z_post[idx,1]*sum.match$SE_tar[idx]
    sum.match$BETA = post_beta_tar
    prs.file = sum.match %>% 
      filter(CHR==j) %>% 
      rename(A1 = A1_tar) %>% 
      select(SNP,A1,BETA)
    prs.eb= prs.file
    
    BETA_eur = prs.eur$BETA
    BETA_tar = prs.tar$BETA
    BETA_tar[is.na(BETA_tar)] = 0
    BETA_eb = prs.eb$BETA
    BETA_eb[is.na(BETA_eb)] = 0
    prs.file = cbind(prs.eur,BETA_tar,BETA_eb)
    
    
    write.table(prs.file,file = paste0(temp.dir,"best_eur_coeff_chr_",j),col.names = T,row.names = F,quote=F)
    
    
    res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink2 --score-col-nums 3,4,5 --threads 2 --score ",temp.dir,"best_eur_coeff_chr_",j," header no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"best_eur_prs_chr_",j))
    
  #}

}






# result.matrix <- matrix(0,3,1)
# for(l in 1:3){
#   #for(i1 in 1:2){
#     sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep)))
#     idx <- which(sum.data$P<=5E-08)
# 
#     result.matrix[l,i1] <- length(idx)
#   #}
# }
