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
j = as.numeric(args[[3]])
l = as.numeric(args[[4]])
m = as.numeric(args[[5]])
i1 = as.numeric(args[[6]])

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

system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.* ",temp.dir))
system(paste0("ls ",temp.dir))

library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/")

pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n.snp.mat <- matrix(0,length(pthres),4)


#for(l in 1:3){
print(l)
# for(m in 1:4){

print(m)
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>% 
  mutate(sd_eur=beta_eur/STAT)
  select(SNP,A1,beta_eur,sd_eur,peur) %>% 
  rename(A1.EUR = A1)
r2_vec = c(0.01,0.05,0.1,0.2,0.5)
wc_base_vec = c(50,100,200,500)
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com <- left_join(sum.data,summary.eur.select,by="SNP")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id,all_of(eth))
summary.com.update = left_join(summary.com,snp.infor,by="SNP")
mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(mega.list) = "rs_id"
summary.com.match = inner_join(summary.com.update,mega.list,by="rs_id")

idx <- which(summary.com.match$A1!=summary.com.match$A1.EUR)
summary.com.match$A1.EUR[idx] <- summary.com.match$A1[idx]
summary.com.match$beta_eur[idx] <- -summary.com.match$beta_eur[idx]

summary.com.match = summary.com.match %>% 
  rename(MAF = eth[i]) %>% 
  mutate(sd_tar=BETA/STAT,
         beta_st = BETA*sqrt(2*MAF*(1-MAF)),
         sd_st = sd_tar*sqrt(2*MAF*(1-MAF)),
         beta_eur_st = beta_eur*sqrt(2*EUR*(1-EUR)),
         sd_eur_st = sd_eur*sqrt(2*EUR*(1-EUR)))
#estimate the prior
load(paste0(out.dir,eth[i],"/r2.list_rho_two_way_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))

p.k1 =r2.list[[4]][[3]]
p.k2 = r2.list[[4]][[4]]
r2_ind = r2.list[[4]][[5]]
w_ind = r2.list[[4]][[6]]
LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
clump.snp <- LD[,3,drop=F] 
summary.com.prior = left_join(clump.snp,summary.com.match,by="SNP") %>% 
  filter(peur<p.k1|
           P<p.k2)
prior.sigma = cov(cbind(summary.com.prior$beta_st,
                        summary.com.prior$beta_eur_st),use="complete.obs")

beta_st = summary.com.match$beta_st
var_st = summary.com.match$sd_st^2
MAF = summary.com.match$MAF
beta_eur_st = summary.com.match$beta_eur_st
var_eur_st = summary.com.match$sd_eur_st^2
MAF_EUR = summary.com.match$EUR

post_beta_tar = beta_st
PostBeta <- function(beta,Sigma,Sigma0,MAF.train,
                     MAF.ref){
  n <- length(beta)
  beta_post <- solve(solve(Sigma)+solve(Sigma0))%*%(solve(Sigma)%*%beta)
  beta_post[1] <- beta_post[1]/sqrt(2*MAF.train*(1-MAF.train))
  beta_post[2] <- beta_post[2]/sqrt(2*MAF.ref*(1-MAF.ref))
  return(beta_post)
}


for(m_i in 1:nrow(summary.com.match)){
  #if(m_i%%1000)print(m_i)
  Sigma = diag(c(var_st[m_i],
                 var_eur_st[m_i]))
  beta = c(beta_st[m_i],
           beta_eur_st[m_i])
  MAF.train = MAF[m_i]
  MAF.ref = MAF_EUR[m_i]
  
  if(is.na(det(Sigma))==0){
    if(det(Sigma)!=0){
      post_beta_tar[m_i] = PostBeta(beta,Sigma,prior.sigma,MAF.train,MAF.ref)[1]   
    }
    
  }
}



summary.com.match$BETA = post_beta_tar
summary.com  = summary.com.match
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    
    
    #read LD clumped SNPs
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD[,3,drop=F] 
    #read the target ethnic group summary level statistics
    #combine the statistics with SNPs after clumping
    prs.all <- left_join(clump.snp,summary.com,by="SNP") 
    colSums(is.na(prs.all))
    for(k1 in 1:length(pthres)){
      #keep al the SNPs with peur pass the threshold
      prs.all.temp = prs.all
      idx <- which(prs.all.temp$peur<=pthres[k1])
      prs.all.temp$P[idx] = 1E-20
      prs.file <- prs.all.temp %>% filter(CHR==j) %>% 
        select(SNP,A1,BETA,P)
      colSums(is.na(prs.file))
      write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
      
      p.value.file <- prs.all.temp %>% filter(CHR==j) %>% 
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
      res = system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"chr",j,".mega --exclude ",old.out.dir,eth[i],"/duplicated.id  --out ",temp.dir.prs,"prs_eb_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
      print("step2 finished")
      #system(paste0("ls ",temp.dir.prs))
      if(res==2){
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








