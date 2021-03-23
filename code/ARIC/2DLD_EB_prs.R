#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is trait
#j is the number of chromsome
args = commandArgs(trailingOnly = T)

j = as.numeric(args[[1]])
#l = as.numeric(args[[3]])
#m = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
#i_rep = 1

#j = as.numeric(args[[3]])
library(dplyr)
library(data.table)

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
n_pthres  = length(pthres)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
i = 2
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)


#generate q_range file is a one time job
# for(l in 1:3){
#   temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
#   q_range = data.frame(
#     paste0("p_value_",c(1:n_pthres)),
#     rep(0,n_pthres),
#     pthres,stringsAsFactors = F)
#   write.table(q_range,file = paste0(temp.dir,"2Dq_range_file"),row.names = F,col.names = F,quote=F)
# 
# }


for(l in 1:3){
  setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
  data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
  temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
  out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
  out.dir.eur = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[1],"/")
  summary.eur = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[1],"/sumdata/training-GWAS-formatted.txt")))
  summary.eur.select = summary.eur %>% 
    mutate(p_eur = PVAL,
           beta_eur = BETA,
           SNP = SNP_ID,
           se_eur = SE,
           A1_eur = REF) %>% 
    select(SNP,beta_eur,se_eur,A1_eur,p_eur)
  # colnames(summary.eur)[9] = "peur"
  # colnames(summary.eur)[7] = "beta_eur"
  # summary.eur.select = summary.eur %>% 
  #
  
  #load target gwas summary statistics
  sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
  sum.data = sum.data %>% 
    mutate(BP=POS,
           SNP = SNP_ID,
           beta_tar = BETA,
           A1_tar = REF,
           p_tar = PVAL,
           se_tar = SE) %>% 
    select(SNP,beta_tar,se_tar,A1_tar,p_tar)
  
  # filter.data <- sum.data %>% 
  #   filter(SNP_ID=="rs1967017")
  # 
  sum.com <- left_join(sum.data,summary.eur.select,by="SNP")
  colSums(is.na(sum.com))
  #align alleles
  sum.com = sum.com %>% 
    mutate(beta_eur=ifelse(A1_tar==A1_eur,
                           beta_eur,
                           -beta_eur),
           A1_eur =ifelse(A1_eur==A1_tar,A1_eur,A1_tar)) %>% 
    mutate(z_stat_eur = beta_eur/se_eur,
           z_stat_tar = beta_tar/se_tar)
  EstimatePrior <- function(beta_tar,sd_tar,
                            beta_eur,sd_eur){
    beta_tar = as.numeric(beta_tar)
    sd_tar = as.numeric(sd_tar)
    beta_eur = as.numeric(beta_eur)
    sd_eur = as.numeric(sd_eur)
    z_tar = beta_tar/sd_tar
    z_eur = beta_eur/sd_eur
    z_mat <-na.omit(cbind(z_tar,z_eur))
    
    
    prior.mat <- cov(z_mat)-diag(2)
    return(prior.mat)
  }
  EBpost <- function(beta_tar,sd_tar,
                     beta_eur,sd_eur,EBprior){
    
    beta_tar = as.numeric(beta_tar)
    sd_tar = as.numeric(sd_tar)
    beta_eur = as.numeric(beta_eur)
    sd_eur = as.numeric(sd_eur)
    prior.sigma = EBprior
    z_tar = beta_tar/sd_tar
    z_eur = beta_eur/sd_eur
    z_mat <-as.matrix(cbind(z_tar,z_eur))
    sd_mat =  as.matrix(cbind(sd_tar,sd_eur))
    post.sigma = solve(solve(prior.sigma)+diag(2))
    
    z_mat_post = z_mat
    
    p <- ncol(z_mat)
    
    for(k in 1:nrow(z_mat)){
      if(k%%10000==0){print(paste0(k," SNPs completed"))}
      z_temp = z_mat[k,]
      
      #find out nonmissing component
      
      idx <- which(!is.na(z_temp))
      if(length(idx)<p){
        z_temp <- z_temp[idx]
        
        post.sigma_temp = post.sigma[idx,idx,drop=F]
        z_post = post.sigma_temp%*%z_temp
      }else{
        z_post =post.sigma%*%z_temp
      }   
      
      z_mat_post[k,idx] = z_post
    }
    beta_mat_post = z_mat_post
    beta_mat_post[,1] =z_mat_post[,1]*sd_tar
    beta_mat_post[,2] =z_mat_post[,2]*sd_eur
    return(beta_mat_post)
  }
  
  
  #load the best 2DLD cutoff to estimate prior
  load(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.2DLD.rdata"))
  r2.result <- ARIC.result.2DLD[[2]]
  
  r2.result.filter <- r2.result %>% 
    filter(trait.vec==trait[l]) %>% 
    filter(rer2.vec.test.prs==max(rer2.vec.test.prs))
  p.k1 =   r2.result.filter$p.eur.vec
  p.k2 = r2.result.filter$p.tar.vec
  r_ind = r2.result.filter$r2.vec
  w_ind = r2.result.filter$wc.vec
  load(paste0(temp.dir,"/2DLD_LD_clump_rind_",r_ind,"_wcind_",w_ind,".rdata"))
  summary.com.prior = left_join(LD,sum.com,by="SNP") %>% 
    filter(p_eur<p.k1|
             p_tar<p.k2) %>% 
    mutate(z_stat_eur = beta_eur/se_eur,
           z_stat_tar = beta_tar/se_tar)
  
  beta_tar <- summary.com.prior$beta_tar
  sd_tar <- summary.com.prior$sd_tar
  beta_eur <- summary.com.prior$beta_eur
  sd_eur <- summary.com.prior$sd_eur
  
  EBprior = EstimatePrior(beta_tar,sd_tar,
                          beta_eur,sd_eur)
  
  
  
  
  
  
  prior.sigma = cov(cbind(summary.com.prior$z_stat_tar,
                          summary.com.prior$z_stat_eur),use="complete.obs")-diag(2)
  
  post.sigma = solve(solve(prior.sigma)+diag(2))
  z_stat_tar = sum.com$z_stat_tar
  z_stat_eur = sum.com$z_stat_eur
  z_mat = cbind(z_stat_tar,z_stat_eur)
  z_post = z_mat%*%post.sigma
  
  sum.com = sum.com %>% 
    mutate(z_post_tar = z_post[,1]) %>% 
    mutate(post_beta_tar = ifelse(is.na(z_post_tar),beta_tar,z_post_tar*se_tar))
    
  
  #load bim file to match all the SNPs with different strand
  #in summary data and ARIC data
  bim <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j,".bim")))
  colnames(bim)[2] <- "SNP"
  sum.com.match <- left_join(bim,sum.com,by="SNP") %>% 
    mutate(A1 = A1_tar,
           BETA = post_beta_tar,
           P = p_tar)
  #A1 is in summary data
  #V5 and V6 is from genotype data
  #if A1 is not either V5 or V6, flip strand
  sum.com.match = sum.com.match %>% 
    mutate(A1 = case_when(
      (A1!=V5&A1!=V6)&A1=="A" ~"T",
      (A1!=V5&A1!=V6)&A1=="C" ~"G",
      (A1!=V5&A1!=V6)&A1=="T" ~"A",
      (A1!=V5&A1!=V6)&A1=="G" ~"C",
      TRUE ~ A1
    ))
  
  # idx <- which((sum.com.match$A1_new!=sum.com.match$V5)&
  #                (sum.com.match$A1_new!=sum.com.match$V6))
  
  for(r_ind in 1:length(r2_vec)){
    wc_vec = wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      
      LD <- as.data.frame(fread(paste0(temp.dir,"2DLD_clump_chr_",j,"_",r_ind,"_",w_ind,".clumped")))
      clump.snp <- LD[,2,drop=F]
      
      prs.all <- left_join(clump.snp,sum.com.match,by="SNP")
      colSums(is.na(prs.all))
      
      
      
      
      
      
      for(k1 in 1:length(pthres)){
        
        #keep al the SNPs with peur pass the threshold
        prs.file <- prs.all %>% 
          mutate(P = replace(P,p_eur<=pthres[k1],1E-20))%>%
          select(SNP,A1,BETA,P)
        write.table(prs.file,file = paste0(temp.dir,"EBprs_coeff_chr_",j),col.names = T,row.names = F,quote=F)
        
        p.value.file = prs.file %>% 
          select(SNP,P)
        
        write.table(p.value.file,file = paste0(temp.dir,"EBp_value_chr_",j),col.names = T,row.names = F,quote=F)
        if(nrow(prs.file)>0){
          res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --q-score-range ",temp.dir,"2Dq_range_file ",temp.dir,"EBp_value_chr_",j," header --threads 2 --score ",temp.dir,"EBprs_coeff_chr_",j," header no-sum no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"EBprs_chr_",j,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
          #res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_chr_",j," header --threads 2 --score ",temp.dir,"prs_coeff_chr_",j," header no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"prs_chr_",j))
          print("step2 finished")
          #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
          if(res==2){
            stop()
          }
          
          
        }
      }
    }
  }
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
