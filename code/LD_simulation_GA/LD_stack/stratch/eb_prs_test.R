#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
i_rep = as.numeric(args[[1]])
i = as.numeric(args[[2]])
l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
#j = as.numeric(args[[3]])
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')

system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))
system(paste0("ls ",temp.dir))

library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/")

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#n.snp.mat <- matrix(0,length(pthres),4)


#for(l in 1:3){
print(l)
# for(m in 1:4){

print(m)
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>% 
  mutate(sd_eur=beta_eur/STAT) %>% 
  select(SNP,A1,beta_eur,sd_eur,peur) %>% 
  rename(A1.EUR = A1)
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com.match <- left_join(sum.data,summary.eur.select,by="SNP")

idx <- which(summary.com.match$A1!=summary.com.match$A1.EUR)
summary.com.match$A1.EUR[idx] <- summary.com.match$A1[idx]
summary.com.match$beta_eur[idx] <- -summary.com.match$beta_eur[idx]

summary.com.match = summary.com.match %>%
  mutate(beta_tar = BETA,
         sd_tar=BETA/STAT,
         z_stat_tar = STAT,
         z_stat_eur = beta_eur/sd_eur)
#estimate the prior
load(paste0(out.dir,eth[i],"/r2.list_rho_2DLD_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))

p.k1 =r2.list[[3]][[3]]
p.k2 = r2.list[[3]][[4]]
r_ind = r2.list[[3]][[5]]
w_ind = r2.list[[3]][[6]]
LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
#LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
clump.snp <- LD[,1,drop=F] 

EstimatePrior <- function(beta_tar,sd_tar,
                          beta_eur,sd_eur){
  beta_mat <-na.omit(cbind(beta_tar,beta_eur))
  var_mat <- na.omit(cbind(sd_tar^2,sd_eur^2))
  
  prior.mat <- (crossprod(beta_mat)-diag(colSums(var_mat)))/(nrow(beta_mat)-1)
  return(prior.mat)
}


EBpost <- function(beta_tar,sd_tar,
                   beta_eur,sd_eur,EBprior){
  
 
  prior.mat.inverse = solve(EBprior)
  beta_mat <-cbind(beta_tar,beta_eur)
  var_mat_inverse <- cbind(1/sd_tar^2,1/sd_eur^2)
  
  beta_mat_post = beta_mat
  
  p <- ncol(beta_mat)
  
  for(k in 1:nrow(beta_mat)){
    if(k%%10000==0){print(paste0(k," SNPs completed"))}
    beta_temp = beta_mat[k,]
    var_mat_inverse_temp = var_mat_inverse[k,]
    #find out nonmissing component
    
    idx <- which(!is.na(beta_temp))
      if(length(idx)<p&length(idx)>0){
      beta_temp <- beta_temp[idx]
      var_mat_inverse_temp = var_mat_inverse_temp[idx,drop=F]
      prior.mat.inverse_temp = prior.mat.inverse[idx,idx,drop=F]
      beta_post = solve(var_mat_inverse_temp+prior.mat.inverse_temp)%*%beta_temp
    }else{
      beta_post = solve(var_mat_inverse_temp+prior.mat.inverse)%*%beta_temp
    }   
    beta_mat_post[k,idx] = beta_post
    
  }
  return(beta_mat_post)
}



summary.com.prior = left_join(clump.snp,summary.com.match,by="SNP") %>% 
  filter(peur<p.k1|
           P<p.k2)
beta_tar <- summary.com.prior$beta_tar
sd_tar <- summary.com.prior$sd_tar
beta_eur <- summary.com.prior$beta_eur
sd_eur <- summary.com.prior$sd_eur

EBprior = EstimatePrior(beta_tar,sd_tar,
                 beta_eur,sd_eur)

beta_tar <- summary.com.match$beta_tar
sd_tar <- summary.com.match$sd_tar
beta_eur <- summary.com.match$beta_eur
sd_eur <- summary.com.match$sd_eur

post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)

post_beta_tar = post_beta_mat[,1,drop=F]



summary.com.match$BETA = post_beta_tar

summary.com  = summary.com.match

#remove duplicated snp
idx <- which(duplicated(summary.com$SNP))
if(length(idx)!=0){
  summary.com = summary.com[-idx,]
}
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    
    
    #read LD clumped SNPs
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD[,1,drop=F] 
    #read the target ethnic group summary level statistics
    #combine the statistics with SNPs after clumping
    prs.all <- left_join(clump.snp,summary.com,by="SNP") 
    #idx <- which(summary.com$SNP=="rs201335322:154415264:C:T")
    colSums(is.na(prs.all))
  
  
    for(k1 in 1:length(pthres)){
      #keep al the SNPs with peur pass the threshold
      prs.all.temp = prs.all
      idx <- which(prs.all.temp$peur<=pthres[k1])
      prs.all.temp$P[idx] = 1E-20
      prs.file <- prs.all.temp %>% 
        select(SNP,A1,BETA,P)
      colSums(is.na(prs.file))
      write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
      
      #idx <- which(prs.file$SNP=="rs4806716:54639868:G:A")
      
      p.value.file <- prs.all.temp %>% 
        select(SNP,P)
      write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
      n_pthres = length(pthres)
      q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
      temp = 1
      
      for(k2 in 1:length(pthres)){
        prs.file.sub <- prs.file %>% 
          filter(P<=pthres[k2]) %>% 
          select(SNP,A1,BETA,P)
        if(nrow(prs.file.sub)>0){
          q_range[temp,1] = paste0("p_value_",k2)
          q_range[temp,3] = pthres[k2]
          temp = temp+1
        }
      }
      q_range = q_range[1:(temp-1),]
      write.table(q_range,file = paste0(temp.dir.prs,"q_range_file"),row.names = F,col.names = F,quote=F)
      res = system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"all_chr_test.mega --exclude ",old.out.dir,eth[i],"/duplicated.id  --out ",temp.dir.prs,"prs_eb_test_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
      print("step2 finished")
      #system(paste0("ls ",temp.dir.prs))
      if(res==3){
        stop()
      }
      res = system(paste0("mv ",temp.dir.prs,"*.profile ",out.dir,eth[i],"/prs/"))
      if(res==2){
        stop()
      }
      system(paste0("rm -rf ",temp.dir.prs))
      dir.create(paste0(temp.dir.prs),showWarnings = FALSE)
      
      
      #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
    }
    
  }
  
}


#pthres <- c(1E-10,1E-09,5E-08,1E-07,2.5E-07,5E-07,7.5E-07,1E-06,2.5E-06,5E-06,7.5E-06,1E-05,2.5e-05,5E-05,7.5e-05,1E-04,2.5E-04,5E-04,7.5E-04,1E-03)

system(paste0("rm -rf ", temp.dir))

#}








