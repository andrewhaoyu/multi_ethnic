#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
i_rep = 10
i =2
l = 1
m = 1
i1 = 1
i_c = 1
k1 = as.numeric(args[[1]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/test_chip/"
#j = as.numeric(args[[3]])
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
if(i_c==2){
  #hapmap3
  system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))
  system(paste0("ls ",temp.dir))
  
}else{
  system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.tag.* ",temp.dir))
}
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
load(paste0(out.dir,eth[i],"/r2.list_rho_2DLD_i_c_",i_c,"_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
idx <- which.max(result.data[,1])
name = result.data$name[idx]
name_split = strsplit(name,"_")
p.k1 = pthres[as.numeric(name_split[[1]][[8]])]
p.k2 = pthres[as.numeric(name_split[[1]][[10]])]
r_ind = as.numeric(name_split[[1]][[3]])
w_ind = as.numeric(name_split[[1]][[6]])
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

LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_i_c_",i_c,"rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))

summary.com.prior = left_join(LD,summary.com.match,by="SNP") %>% 
  filter(peur<p.k1|
           P<p.k2)
beta_tar <- summary.com.prior$beta_tar
sd_tar <- summary.com.prior$sd_tar
beta_eur <- summary.com.prior$beta_eur
sd_eur <- summary.com.prior$sd_eur

EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)





rs.id.list = list()
temp = 1
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_i_c_",i_c,"rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    rs.id.list[[temp]] = LD   
    temp = temp + 1
  }
}
unique.id = unique(rbindlist(rs.id.list))

unique.infor = left_join(unique.id,summary.com.match,by="SNP")
idx <- which(duplicated(unique.infor$SNP))
if(length(idx)!=0){
  unique.infor = unique.infor[-idx,]
}





beta_tar <- unique.infor$beta_tar
sd_tar <- unique.infor$sd_tar
beta_eur <- unique.infor$beta_eur
sd_eur <- unique.infor$sd_eur
post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)
post_beta_mat[is.na(post_beta_mat)] = 0
n.col = length(r2_vec)*length(wc_vec)
beta_mat = matrix(rep(post_beta_mat,n.col),ncol =n.col*ncol(post_beta_mat))
temp = 0
names = rep("c",ncol(beta_mat))
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    #summary.com  = summary.com.match
    
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_i_c_",i_c,"rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    idx <- which(unique.id$SNP%in%LD$SNP==F)
    beta_mat[idx,(1:ncol(post_beta_mat))+temp] = 0
    names[(1:ncol(post_beta_mat))+temp] = paste0(colnames(post_beta_mat),"_r_ind_",r_ind,"_w_ind_",w_ind)
    temp = temp + ncol(post_beta_mat)
    
  }
}

colnames(beta_mat) = names
prs.file = data.frame(SNP = unique.id,A1 = unique.infor$A1,beta_mat)
n.col = ncol(prs.file)
write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)

#create q-range file
n_pthres = length(pthres)
q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
temp = 1

for(k2 in 1:length(pthres)){
  
  q_range[temp,1] = paste0("p_value_",k2)
  q_range[temp,3] = pthres[k2]
  temp = temp+1
}
q_range = q_range[1:(temp-1),]
write.table(q_range,file = paste0(temp.dir.prs,"q_range_file"),row.names = F,col.names = F,quote=F)

#create p-value file
p.value.file = data.frame(SNP = unique.id,P = unique.infor$P)
#for(k1 in 1:length(pthres)){
  #keep al the SNPs with peur pass the threshold
  idx <- which(unique.infor$peur<=pthres[k1])
  p.value.file$P[idx] = 1E-20
  
  write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
  
  
  if(i_c==2){
    res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                        "--q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header ",
                        "--threads 2 ",
                        "--score-col-nums 3-",n.col," ",
                        "--score ",temp.dir.prs,"prs_file header no-mean-imputation ",
                        "--bfile ",temp.dir,"all_chr_test.mega ",
                        "--exclude ",old.out.dir,eth[i],"/duplicated.id  ",
                        "--out ",temp.dir.prs,"prs_eb_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"p_value_",k1))  
    
  }else{
    res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                        "--q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header ",
                        "--threads 2 ",
                        "--score-col-nums 3-",n.col," ",
                        "--score ",temp.dir.prs,"prs_file header no-mean-imputation ",
                        "--bfile ",temp.dir,"all_chr_test.tag ",
                        "--exclude ",old.out.dir,eth[i],"/duplicated.id  ",
                        "--out ",temp.dir.prs,"prs_eb_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"p_value_",k1))  
    
  }
  
  print("step2 finished")
  #dup <- fread(paste0(old.out.dir,eth[i],"/duplicated.id"),header = F)
  
  
  #system(paste0("ls ",temp.dir.prs))
  if(res==3){
    stop()
  }
  #res = system(paste0("mv ",temp.dir.prs,"*.profile ",out.dir,eth[i],"/prs/"))
  if(res==2){
    stop()
  }
  
  #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
  
  
#}

#combine all the prs
prs.list = list()
temp = 1
#for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    prs.temp = fread(paste0(temp.dir.prs,"prs_eb_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"p_value_",k1,".p_value_",k2,".sscore"))
    prs.list[[temp]] = prs.temp[,5:ncol(prs.temp)]
    colnames(prs.list[[temp]]) = paste0(names,"_","p_value_",k1,"_p_value_",k2)
    temp = temp + 1
  }
#}
prs.mat = cbind(prs.temp[,1:2],bind_cols(prs.list))
write.table(prs.mat,file = paste0(out.dir,eth[i],"/prs/prs_eb_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".profile"))



prs.list = list()
temp = 1
#for(k1 in 1:length(pthres)){
for(k1 in 1:length(pthres)){
  prs.temp = fread(paste0(out.dir,eth[i],"/prs/prs_eb_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".profile"))
  prs.list[[temp]] = prs.temp[,4:ncol(prs.temp)]
#  colnames(prs.list[[temp]]) = paste0(names,"_","p_value_",k1,"_p_value_",k2)
  temp = temp + 1
}
#}
prs.mat = cbind(prs.temp[,2:3],bind_cols(prs.list))
write.table(prs.mat,file = paste0(out.dir,eth[i],"/prs/prs_eb_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".profile"))




