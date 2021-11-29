#merge the prs by chromosome to one for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#get the number of SNPs used in weighted-PRS method
method <- c("Weighted-PRS")
i_rep_saved = i_rep
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
setwd("/data/zhangh24/multi_ethnic/")


total <- 1
eth.vec <- rep("c",total)
n.snp.result <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
ga_vc <- rep(0,total)
method_vec <- rep("c",total)
temp= 1
# for(i in 2:5){
#   for(i1 in 1:5){
#     for(l in 1:3){
      #load the phenotype file
      y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
      y <- y[,2+(1:n.rep)]
      n <- nrow(y)
      y_test_mat <- y[(100000+1):110000,]
      y_vad_mat <- y[(110001):120000,]
      #for(m in 1:4){
        print(m)
        
       
        
        # for(i_rep in 1:n.rep){
          
      
          #get the number of SNPs in EUR PRS
          #load EUR r2 result
          out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
          load(paste0(out.dir,"LD.clump.result.CT.rdata"))
          
          #load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
          #keep the EUR results with sample size at 100,000 (m = 4)
          r2.mat <- as.data.frame(LD.result.list[[2]]) %>% 
            filter(eth_vec=="EUR"&
                     m_vec==4&
                     l_vec==l&
                     ga_vec==i1) %>% 
            filter(r2.vec.test==max(r2.vec.test))
          #get the best performance eur prs p-value threshold
          p.cut = r2.mat$pthres_vec
          r_ind = 1
          w_ind = 1
          LD <- as.data.frame(fread(paste0(out.dir,eth[1],"/LD_clump_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
          clump.snp <- LD[,3,drop=F] 
          sum.data <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
          colnames(sum.data)[2] <- "SNP"
          prs.all <- left_join(clump.snp,sum.data,by="SNP") 
          prs.file <- prs.all %>% filter(P<=p.cut) %>% 
            mutate(SE_eur = BETA/STAT) %>% 
            select(SNP)
          prs.eur = prs.file
          
          #get the number of the target population
          filename <- paste0(out.dir,eth[i],"/prs/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".sscore")
          prs.temp <- as.data.frame(fread(filename))
          #load AUC performance for target
          load(paste0(out.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_GA_",i1))
          best.auc.result <- as.data.frame(r2.list[[2]]) %>% 
            filter(r2.vec.test==max(r2.vec.test))
          r_ind  = best.auc.result$r2_ind_vec
          w_ind = best.auc.result$wc_ind_vec
          k =which(pthres==best.auc.result$pthres_vec)
          sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
          LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
          clump.snp <- LD[,3,drop=F]           
          pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
          prs.file <- left_join(clump.snp,sum.data,by="SNP") %>% 
            filter(P<pthres[k]) %>% 
            select(SNP)
          prs.tar = prs.file
          prs.snp = unique(rbind(prs.eur,prs.tar))
      
          
          
          n.snp.result <- nrow(prs.snp)
          
          
          
        # }
        
        
        eth.vec[temp] <- eth[i]
        l_vec[temp] <- l
        m_vec[temp] <- m
        ga_vc[temp] <- i1
        method_vec[temp] <- method
    
        
#       }
#     }
#   }
#   
#   
# }


weightedprs.result <- data.frame(eth.vec,
                                 n.snp.vec = n.snp.result,
                                 l_vec,
                                 m_vec,
                                 ga_vec=ga_vc,
                                 method_vec = method_vec)


#get the number of SNPs in TDLD-EB and TDLD-SLEB
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#n <- 120000

#for(m in 1:1){
i_rep = i_rep_saved
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

r2.vec.test <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
r2.vec.vad <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec1 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec2 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
r2_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
wc_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
prs.mat <- matrix(0,n.test+n.vad,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
temp = 1
snp.list = list()
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
        snp.list[[temp]] = prs.all %>% select(SNP)
        if(nrow(prs.all)>0){
          filename <- paste0(out.dir,eth[i],"/prs/prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
          
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE
          prs.test <- prs.score[(1):(n.test)]
          prs.vad <- prs.score[(n.test+1):(n.test+n.vad)]
          #model = lm(y~prs.score)
          y.test = y_test_mat[1:n.test,i_rep]
          y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
          model1 <- lm(y.test~prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          model2 <- lm(y.vad~prs.vad)
          r2.vec.test[temp] = summary(model1)$r.square
          r2.vec.vad[temp] = summary(model2)$r.square
          pthres_vec1[temp] = pthres[k1]
          pthres_vec2[temp] = pthres[k2]
          r2_ind_vec[temp] = r_ind
          wc_ind_vec[temp] = w_ind
          prs.mat[,temp] = prs.score
          temp = temp+1
        }else{
          r2.vec.test[temp] = 0
          r2.vec.vad[temp] = 0
          pthres_vec1[temp] = pthres[k1]
          pthres_vec2[temp] = pthres[k2]
          r2_ind_vec[temp] = r_ind
          wc_ind_vec[temp] = w_ind
          prs.mat[,temp] = 0
          temp = temp+1
        }
        
      }
    }
    
  }
}
result.data <- data.frame(r2.vec.test,r2.vec.vad,
                          pthres_vec1,pthres_vec2,
                          r2_ind_vec,
                          wc_ind_vec)
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

prs.idx = setdiff(c(1:(temp-1)),idx[drop])
prs.mat.new = prs.mat %>% 
  select(-all_of(drop))

select.snp.list = list()
temp = 1
for(k in prs.idx){
  select.snp.list[[temp]] = snp.list[[k]]
  temp = temp+1
}
select.snp = unique(rbindlist(select.snp.list))
n.snp.result = nrow(select.snp)
method_vec = "TDLD-SLEB"
tdldsleb.result <- data.frame(eth.vec,
                                 n.snp.vec = n.snp.result,
                                 l_vec,
                                 m_vec,
                                 ga_vec=ga_vc,
                                 method_vec = method_vec)

#TDLD-EB number of SNPs
idx <- which.max(r2.vec.test)
n.snp.result <- nrow(snp.list[[idx]])
method_vec = "TDLD-EB"
eb.result <- data.frame(eth.vec,
                              n.snp.vec = n.snp.result,
                              l_vec,
                              m_vec,
                              ga_vec=ga_vc,
                              method_vec = method_vec)
#PRS-CSx result
setwd(paste0(out.dir.sum,eth[i],"/prscsx/"))
phi = c("1e-02","1e-04","1e-06")
snp.list = list()
temp = 1
for(v in 1:3){
  #load target population posterior
  prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                     "_pst_eff_a1_b0.5_phi",phi[v],".txt"))
  prs_infor = prs %>% 
    select(V2,V4)
  snp.list[[temp]] = prs_infor %>% select(V2)
  temp = temp + 1
}
#load EUR PRS
for(v in 1:3){
  prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi"
                     ,phi[v],".txt"))
  prs_infor = prs %>% 
    select(V2,V4)
  snp.list[[temp]] = prs_infor %>% select(V2)
  temp = temp + 1
}

select.snp = unique(rbindlist(snp.list))
n.snp.result = nrow(select.snp)
prscsx.result <- data.frame(eth.vec,
                        n.snp.vec = n.snp.result,
                        l_vec,
                        m_vec,
                        ga_vec=ga_vc,
                        method_vec = method_vec)

result = rbind(weightedprs.result,
               tdldsleb.result,
               eb.result,
               prscsx.result)
save(result,file = paste0(out.dir,"n_snp_result_eth_",i,"_rho_",l,"_size_",m,"_ga_",i1,".rdata"))
