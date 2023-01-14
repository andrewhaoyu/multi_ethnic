args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

# for(i in 2:3){
#   for(l in 1:2){
    library(dplyr)
    library(data.table)
    eth <- c("EUR","AFR","AMR")
    trait_vec <-c("height","bmi")
    trait = trait_vec[l]
    out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSX/",eth[i],"/",trait,"/")
    phi = c("1e+00","1e-02","1e-04","1e-06")
    
    sid<-Sys.getenv('SLURM_JOB_ID')
    dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
    temp.dir = paste0('/lscratch/',sid,'/test/')
    
    system(paste0("mkdir ",temp.dir,"ukb"))
    geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
    system(paste0("cp ",geno.data,eth[i],"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
    system(paste0("cp ",geno.data,eth[i],"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
    system(paste0("cp ",geno.data,eth[i],"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))
    v = 1
    #find unique snp set
    #load target prs
    load(paste0(out.dir.prs,"sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
    data_tar = data %>% select(V2, V4)
    #load eur prs
    load(paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
    data_eur = data %>% select(V2, V4)
    snp_set = rbind(data_tar, data_eur) %>% 
      distinct() 
    colnames(snp_set) = c("SNP", "A1")
    
    BETA_mat = matrix(0,nrow(snp_set),length(phi)*2)
    temp = 1
    for(v in 1:4){
      load(paste0(out.dir.prs,"sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
      data_tar = data %>% select(V2, V6) %>% 
        rename(SNP = V2)
      snp_set_temp = left_join(snp_set, data_tar) %>% 
        mutate(BETA = ifelse(is.na(V6),0,V6))
      BETA_mat[, temp] = snp_set_temp$BETA
      temp = temp + 1
      
      load(paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
      data_eur = data %>% select(V2, V6) %>% 
        rename(SNP = V2)
      snp_set_temp = left_join(snp_set, data_eur) %>% 
        mutate(BETA = ifelse(is.na(V6),0,V6))
      BETA_mat[, temp] = snp_set_temp$BETA
      temp = temp + 1
      
    }
    
    prs_file = cbind(snp_set,BETA_mat)
    n_col = ncol(prs_file)
    #file_out = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSx/",eth[i],"/",trait,"")
    write.table(prs_file,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
    
    ref_gene_pred = paste0(temp.dir,"ukb/all_chr")
    res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                        "--score-col-nums 3-",n_col," --threads 2 ",
                        "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                        "--bfile ",ref_gene_pred,
                        " --out ",temp.dir,"PRS"))
    
    prs_mat = fread(paste0(temp.dir,"PRS.sscore"))
    prs_score = prs_mat[,5:ncol(prs_mat)]
    colnames(prs_mat)[2] = "id"
    pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
    pheno_tuning = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth[i],"_tuning.txt")))
    pheno_tuning = pheno_tuning[,1:2]
    covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth[i],"_all_data.txt")))
    pheno_tuning <- left_join(pheno_tuning, covar)
    colnames(pheno_tuning) = c('id','y','sex','age',paste0('pc',1:10))
    pheno_tuning_com = pheno_tuning[complete.cases(pheno_tuning$y),]
    pheno_tuning = left_join(pheno_tuning_com,prs_mat,by = "id")
    model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
    y_tun = model.null$residual
    prs_tun = pheno_tuning[,colnames(prs_score)]
    
    pheno_vad = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth[i],"_validation.txt")))
    pheno_vad = pheno_vad[,1:2]
    pheno_vad <- left_join(pheno_vad, covar)
    colnames(pheno_vad) = c('id','y','sex','age',paste0('pc',1:10))
    pheno_vad_com = pheno_vad[complete.cases(pheno_vad$y),]
    pheno_vad = left_join(pheno_vad_com,prs_mat,by = "id")
    model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
    y_vad = model.null$residual
    prs_vad = pheno_vad[,colnames(prs_score)]
    
    
    r2_vec_tun = rep(0, 4)
    coef_mat = matrix(0, 4, 2)
    for(v in 1:4){
      score1 = prs_tun[,2*v-1]
      score2 = prs_tun[,2*v]
      model = lm(y_tun~score1+score2)
      coef_mat[v,] = coefficients(model)[2:3]
      r2_vec_tun[v] = summary(model)$r.square
    }
    
    max_ind = which.max(r2_vec_tun)
    
    score1 = prs_vad[, 2*max_ind-1]
    score2= prs_vad[, 2*max_ind]
    coef = coef_mat[max_ind,]
    
    prs_vad = cbind(score1, score2)%*%coef
    
    model = lm(y_vad~ prs_vad)
    r2 = summary(model)$r.square
    
    data = data.frame(y = y_vad, x = prs_vad)
    R2Boot = function(data,indices){
      boot_data = data[indices, ]
      model = lm(y ~ x, data = boot_data)
      result = summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 = boot(data = data, statistic = R2Boot, R = 10000)
    
    ci_result = boot.ci(boot_r2, type = "bca")
    
    
    
    r2.result = data.frame(eth = eth[i],
                           trait = trait_vec[l],
                           method = "PRS-CSx",
                           r2 = r2,
                           r2_low = ci_result$bca[4],
                           r2_high = ci_result$bca[5]
    )
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/PRSCSX/",eth[i],"/",trait,"/")
    save(r2.result, file = paste0(out.dir, "prscsx.result"))
#     system(paste0("rm -rf ", temp.dir))
#   }
# }
