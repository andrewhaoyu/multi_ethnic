library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n.rep = 10

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
setwd("/data/zhangh24/multi_ethnic/")


total <- 4*3*4*5
#total <- 4*3*4*1
eth.vec <- rep("c",total)
r2.result <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
ga_vc <- rep(0,total)
method_vec <- rep("PolyPred",total)
temp= 1
#phi = c("1e+00","1e-02","1e-04","1e-06")

mis_vec_list = NULL 
mis_temp = 1
for(i in 2:2){
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass")
  files = dir(file_dir, pattern = paste0(".sscore"),
              full.names = T)
  for(i1 in 1:1){
    for(l in 1:1){
      # for(i in 2:2){
      #   for(i1 in 1:1){
      #     for(l in 1:1){
      
      
      #load the phenotype file
      y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
      y <- y[,2+(1:n.rep)]
      n <- nrow(y)
      y_test_mat <- y[(100000+1):nrow(y),]
      #for(m in 1:4){
      for(m in 1:4){
        #for(q in 1:3){
        
        print(m)
        
        r2.vad.rep <- rep(0,n.rep)
        
        for(i_rep in 1:n.rep){
          r2.test.rep <- rep(0,3)
          coef.mat = matrix(NA,4,2)
          #find the best phi
          
          
          #load target prs
          SBayesR_tar = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_PRS_SBayesR_tar.sscore")
          # if(file_name%in%files){
          PRS <- as.data.frame(fread(paste0(SBayesR_tar)))
          sbayes_tar = PRS$SCORE1_SUM
          sbayes_tar_tun = sbayes_tar[1:10000]
          sbayes_tar_vad = sbayes_tar[10001:20000]
          #
          SBayesR_EUR = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[1],"/polypred/rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1,"_PRS_SBayesR_EUR.sscore")
          PRS <- as.data.frame(fread(paste0(SBayesR_EUR)))
          sbayes_eur = PRS$SCORE1_SUM
          sbayes_eur_tun = sbayes_eur[1:10000]
          sbayes_eur_vad = sbayes_eur[10001:20000]
          #
          polyfun_EUR = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[1],"/polypred/rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1,"_PRS_polyfun.sscore")
          PRS <- as.data.frame(fread(paste0(polyfun_EUR)))
          polyfun_eur = PRS$SCORE1_SUM
          polyfun_eur_tun = polyfun_eur[1:10000]
          polyfun_eur_vad = polyfun_eur[10001:20000]
          
          y.tun = y_test_mat[1:10000,i_rep]
          model = lm(y.tun~sbayes_tar_tun+sbayes_eur_tun+polyfun_eur_tun)
          summary(model)
          coef = coefficients(model)[2:4]
          y.vad = y_test_mat[10001:20000,i_rep]
          prs.tar = cbind(sbayes_tar_vad, sbayes_eur_vad, polyfun_eur_vad)%*%coef
          model = lm(y.vad~prs.tar)
          r2.vad.rep[i_rep] <- summary(model)$r.square
          
          summary(model)
          # }else{
          #   mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
          #   mis_temp = mis_temp + 1
          #   r2.vad.rep[i_rep] = NA
          # }
          #for(k in 1:length(pthres)){
          
          #get the number of
          
          
          
          
          
        }
        
        
        eth.vec[temp] <- eth[i]
        r2.result[temp] <- mean(r2.vad.rep, na.rm = T)
        l_vec[temp] <- l
        m_vec[temp] <- m
        ga_vc[temp] <- i1
        #method_vec[temp] <- method[q]
        temp = temp+1
        
      }
      # }
    }
    
    
  }
}
mis_vec = rbindlist(mis_vec_list)
xpass.result <- data.frame(eth.vec,
                           r2.vec = r2.result,
                           l_vec,
                           m_vec,
                           ga_vec=ga_vc,
                           method_vec = method_vec)

save(xpass.result,file = paste0(out.dir,
                                "xpass.result.rdata"))
