#merge the prs by chromosome to one for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
#i_rep = 2
i1 = 1
i_c = as.numeric(args[[4]])
#for(m in 1:4){
  for(i_rep in 1:10){
    library(dplyr)
    library(data.table)
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
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
    out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/test_chip/"
    y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
    y <- y[,2+(1:n.rep),drop=F]
    n <- nrow(y)
    y_test_mat <- y[(100000+1):nrow(y),,drop=F]
    
    r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
    wc_base_vec = c(50,100)
    
    #get optimal performance on tuning dataset
    r2.vec.test <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
    
    temp = 1
    #all prs is written in one file
    prs.all = fread(paste0(out.dir,eth[i],"/prs/prs_2DLD_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",6,"_wcind_",2,".profile"))
    #drop the first three id columns
    prs.all = as.data.frame(prs.all[,-c(1:3),drop=F])
    
    for(k in 1:ncol(prs.all)){
      prs.score <- prs.all[,k]
      prs.test <- prs.score[(1):(n.test)]
      
      #model = lm(y~prs.score)
      y.test = y_test_mat[1:n.test,i_rep]
      
      model1 <- lm(y.test~prs.test)
      r2.vec.test[k] = summary(model1)$r.square
    }
    name = colnames(prs.all)
    result.data <- data.frame(r2.vec.test,
                              name)
    save(result.data,file = paste0(out.dir,eth[i],"/r2.list_rho_2DLD_i_c_",i_c,"_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    
  }
  
#}
