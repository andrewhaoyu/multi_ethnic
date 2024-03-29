#calculate R2 for PRS using both target and EUR posterior
#merge the prs by chromosome to one for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
#i_rep = 2
i1 = as.numeric(args[[5]])

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)

#n <- 120000

#for(m in 1:1){

n.test <- 10000
n.vad <- n.test
n.rep = 10
#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
setwd("/data/zhangh24/multi_ethnic/")
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
y <- y[,2+(1:n.rep),drop=F]
n <- nrow(y)
y_test_mat <- y[(100000+1):nrow(y),,drop=F]

#for(i_rep in 1:n.rep){
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>%
  select(SNP,beta_eur,peur)
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com <- left_join(sum.data,summary.eur.select,by="SNP")

r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
neth = 2
total = neth*length(pthres)^2*length(r2_vec)*length(wc_base_vec)
r2.vec.test <- rep(0,total)
r2.vec.vad <- rep(0,total)
pthres_vec1 <- rep(0,total)
pthres_vec2 <- rep(0,total)
r2_ind_vec <- rep(0,total)
wc_ind_vec <- rep(0,total)
prs.mat <- matrix(0,n.test+n.vad,total)
temp = 1

for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD
    
    #    colnames(sum.data)[2] <- "SNP"
    
    #for(k in 1:length(pthres)){
    
    prs.clump = left_join(clump.snp,summary.com,by="SNP")
    
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        prs.all <- prs.clump %>%
          filter(peur<=pthres[k1]|
                   P<=pthres[k2])
        if(nrow(prs.all)>0){
          filename <- paste0(out.dir,eth[i],"/prs/prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE
          prs.mat[,temp] = prs.score
          temp = temp+1
        }else{
          prs.mat[,temp] = 0
          temp = temp+1
        }
        #load eur prs
        #take out all the population-specific files
        if(nrow(prs.all)>0&(sum(is.na(prs.all$beta_eur))<nrow(prs.all))){
          filename <- paste0(out.dir,eth[i],"/prs/prs_ebeur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE
          prs.mat[,temp] = prs.score
          temp = temp+1
        }else{
          prs.mat[,temp] = 0
          temp = temp+1
        }
      }
    }
    
  }
}
# result.data <- data.frame(r2.vec.test,r2.vec.vad,
#                           pthres_vec1,pthres_vec2,
#                           r2_ind_vec,
#                           wc_ind_vec)
prs.sum = colSums(prs.mat)
idx <- which(prs.sum!=0)
#drop the prs with all 0
prs.mat <- prs.mat[,idx]
#drop the columns with perfect correlation
prs.mat = as.data.frame(prs.mat)
mtx = cor(prs.mat[1:n.test,])
library(caret)
drop = findCorrelation(mtx,cutoff=0.98)
drop = names(prs.mat)[drop]
prs.mat.new = prs.mat %>% 
  select(-all_of(drop))
num_prs = ncol(prs.mat.new)
save(num_prs,file = paste0(out.dir,eth[i],"/num_prs_rho_ebtest_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#save(r2.list,file = paste0(out.dir,eth[i],"/r2.list_rho_eb_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
