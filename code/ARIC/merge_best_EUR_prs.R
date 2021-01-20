
library(dplyr)
library(data.table)


eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
i = 2


  for(l in 1:3){
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    fam.file <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
    prs.score <- matrix(0,nrow(fam.file),3)
    colnames(prs.score) <- c("eurcoef","tarcoef","ebcoef")
        for(j in 1:22){
  
            
            
            filename <- paste0(temp.dir,"best_eur_prs_chr_",j,".sscore")
            
            prs.temp <- fread(filename)  
            prs.score <- prs.temp[,5:7]* as.numeric(prs.temp[1,3]/2)+prs.score
        }
    write.table(prs.score,file = paste0(out.dir,"/prs_best_eur_prs.profile"),row.names = F,col.names = F,quote=F)
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
