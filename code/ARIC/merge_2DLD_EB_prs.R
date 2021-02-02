
library(dplyr)
library(data.table)

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
#merge clump data
i = 2
r2_vec = c(0.01,0.05,0.1,0.2,0.5)
wc_base_vec = c(50,100)


i = 2
for(l in 1:3){
  temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
  data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
  out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
  files = dir(path = temp.dir,pattern=paste0("EBprs_chr_(.*).profile"),full.names = T)
  for(r_ind in 1:length(r2_vec)){
    for(w_ind in 1:length(wc_base_vec)){
      for(k1 in 1:length(pthres)){
        for(k2 in 1:length(pthres)){
          print(c(k1,k2))
          
          fam.file <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
          prs.score <- rep(0,nrow(fam.file))
          
          for(j in 1:22){
            
            filename <- paste0(temp.dir,"/EBprs_chr_",j,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
            #filename <- paste0(temp.dir,"prs_chr_",j,".p_value_",k,".sscore")
            if(filename%in%files){
              prs.temp <- fread(filename)  
              prs.score <- prs.temp$SCORE*prs.temp$CNT+prs.score
              
            }
            #n.snp.total = n.snp.total+prs.temp$CNT/2
            
          }
          write.table(prs.score,file = paste0(temp.dir,"/EBprs_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,"_",k2,".profile"),row.names = F,col.names = F,quote=F)
          
          
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
