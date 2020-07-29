#Goal: select the causal SNPs
#the heritability are attributable to different proportion of SNPs
#the total heritability are the same
#different proportion of causal SNPs

#load the SNPs information
load("/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
#only select SNPs that don't have duplicated position number
#It will make the extract procedure easier
#the CHR2 EUR file after 700000 are corrupted
#causal SNPs shouldn't come from there
# idx <- which(snp.infor$EUR>=0.01&
#                snp.infor$EUR<=0.99&
#                snp.infor$CHR==2&
#                snp.infor$position>=239554914)
# snp.infor <- snp.infor[-idx,]
# save(snp.infor,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")

EUR <- which(snp.infor$EUR>=0.01&
               snp.infor$EUR<=0.99)
AFR <-  which(snp.infor$AFR>=0.01&
                snp.infor$AFR<=0.99)
AMR <-  which(snp.infor$AMR>=0.01&
                snp.infor$AMR<=0.99)
EAS <-  which(snp.infor$EAS>=0.01&
                snp.infor$EAS<=0.99)
SAS <-  which(snp.infor$SAS>=0.01&
                snp.infor$SAS<=0.99)

library(data.table)
EUR.bi <- snp.infor$EUR>=0.01&
               snp.infor$EUR<=0.99
AFR.bi <-  snp.infor$AFR>=0.01&
                snp.infor$AFR<=0.99
AMR.bi <-  snp.infor$AMR>=0.01&
                snp.infor$AMR<=0.99
EAS.bi <-  snp.infor$EAS>=0.01&
                snp.infor$EAS<=0.99
SAS.bi <-  snp.infor$SAS>=0.01&
                snp.infor$SAS<=0.99
eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
#rho is the proportion of causal SNPs among all the SNPs
rho <- c(0.01,0.001,0.0005)
cau.snp.infor.list <- list()

  set.seed(666)
  all.id <- c(1:nrow(snp.infor))
  
  eth.id <- c(1:5)
  #four sub list for four different proportions
  all.cau.list <- list(list(),
                       list(),
                       list())
  temp <- 1
  for(k in 1:5){
    combn.mat <- combn(5,k)
    
    
    for(l in 1:ncol(combn.mat)){
      #select the ethnic groups in this category
      idx.in <- which(eth.id%in%combn.mat[,l])
      #select the ethnic groups not in this category
      idx.out <- which(eth.id%in%combn.mat[,l]==F)
     
          if(length(idx.out)>=1){
            #sample the smaller proportion from the larger proportion
            #this way i don't need to reextract the SNPs
            #i just need to extract the 10% causal SNPs
            #then the 1%... are all subsets of the 10%
            i = 1
            idx <- which(rowSums(eth.bi[,idx.in,drop=F])==length(idx.in)&
                           rowSums(eth.bi[,idx.out,drop=F])==0)
            cau.idx <- sample(idx,round(length(idx)*rho[i],0),replace=F)
            all.cau.list[[i]][[temp]] <- as.data.frame(cau.idx)
            for(i in 2:length(rho)){
              cau.idx.new <- sample(cau.idx,round(length(cau.idx)*rho[i]/rho[1],0),replace=F)
              all.cau.list[[i]][[temp]] <- as.data.frame(cau.idx.new)
            }
            
  
            temp  = temp + 1
          }else{
            i = 1
            idx <- which(rowSums(eth.bi[,idx.in,drop=F])==length(idx.in))
            
            cau.idx <- sample(idx,round(length(idx)*rho[i],0),replace=F)
            all.cau.list[[i]][[temp]] <- as.data.frame(cau.idx)
            for(i in 2:length(rho)){
              cau.idx.new <- sample(cau.idx,round(length(cau.idx)*rho[i]/rho[1],0),replace=F)
              all.cau.list[[i]][[temp]] <- as.data.frame(cau.idx.new)
            }
            
            temp  = temp + 1 
          }
          
        }
      
      
    }
  
  cau.snp.infor.list <- list()
  for(i in 1:length(rho)){
    cau.id <- rbindlist(all.cau.list[[i]])
    cau.idx <- as.numeric(as.character(as.data.frame(cau.id)[,1]))
    cau.snp.infor.list[[i]] <- snp.infor[cau.idx,]
    
  }

  
  
  
  
  
  save(cau.snp.infor.list,file = "/data/zhangh24/multi_ethnic/result/LD_simulation_new/cau.snp.infor.list.rdata")
  



