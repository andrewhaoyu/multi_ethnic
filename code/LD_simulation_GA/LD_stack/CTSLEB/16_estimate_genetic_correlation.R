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

LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
clump.snp <- LD
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
beta_tar = as.numeric(beta_tar)
sd_tar = as.numeric(sd_tar)
beta_eur = as.numeric(beta_eur)
sd_eur = as.numeric(sd_eur)
prior.sigma = EstimatePrior(beta_tar,sd_tar,beta_eur,sd_eur)
gr = prior.sigma[1,2]/(prior.sigma[1,2]*prior.sigma[2,2])
out.dir = "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"

save(gr, file = paste0(out.dir,eth[i],"/gr_esitmate_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#}








