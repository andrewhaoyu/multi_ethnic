#merge the prs by chromosome to one for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#m = as.numeric(args[[3]])
i_rep = as.numeric(args[[3]])
#i_rep = 2
i1 = as.numeric(args[[4]])

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000

for(m in 1:1){
  
  n.test <- 10000
  n.vad <- n.test
  n.rep = 3
  #r2 mat represent the r2 matrix for the testing dataset
  #column represent the ethnic groups
  #row represent different p-value threshold
  cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  setwd("/data/zhangh24/multi_ethnic/")
  out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
  
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
  
  r2_vec = c(0.01,0.05,0.1,0.2,0.5)
  wc_base_vec = c(50,100,200,500)
  for(r_ind in 1:length(r2_vec)){
    wc_vec = wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      print(c(r_ind,w_ind))
      
      LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
      clump.snp <- LD[,3,drop=F]  
      
      colnames(sum.data)[2] <- "SNP"
      
      #for(k in 1:length(pthres)){
      
      prs.clump = left_join(clump.snp,summary.com,by="SNP")
      
      for(k1 in 1:length(pthres)){
        for(k2 in 1:length(pthres)){
        prs.all <- prs.clump %>% 
          filter(peur<=pthres[k1]|
                   P<=pthres[k2])
        if(nrow(prs.all)>0){
          n.snp.total <- 0
          prs.score <- rep(0,n.test+n.vad)
          for(j in 1:22){
            
            #get the number of
            idx <- which(prs.all$CHR==j)
            n.snp.total = n.snp.total+length(idx)
            if(length(idx)>0){
              
              
              filename <- paste0(out.dir,eth[i],"/prs/prs_chr_",j,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".p_value_",k,".profile")
              
              prs.temp <- fread(filename)  
              prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
            }
          }
          write.table(prs.score,file = paste0(out.dir,eth[i],"/prs/prs_pvalue_",k,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,"p_value,",k2,".profile"),row.names = F,col.names = F,quote=F)
        }
        }
      }
    }
  }
}
#}
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
