
library(dplyr)
library(data.table)

pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
#merge clump data

for(i in 1:2){
  for(l in 1:3){
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    
    LD.list = list()
    for(j in 1:22){
      
      LD <- as.data.frame(fread(paste0(out.dir,"/LD_clump_chr_",j,".clumped")))
      clump.snp <- LD[,3,drop=F]  
      LD.list[[j]] = clump.snp
    }
    LD = rbindlist(LD.list)
    write.table(LD,file = paste0(out.dir,"/LD_clump.clumped"),row.names = F,col.names = T)
  }}

for(i in 1:2){
  for(l in 1:3){
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    
    sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    sum.data = sum.data %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL) 
    # write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
    #             ,col.names = T,row.names = F,quote=F) 
    #prepare association file for plink
    clump.snp <- as.data.frame(fread(paste0(out.dir,"/LD_clump.clumped")))
    prs.clump <- left_join(clump.snp,sum.data,by="SNP")
    
    for(k in 1:length(pthres)){
      prs.all <- prs.clump %>% 
        filter(P<=pthres[k])
      if(nrow(prs.all)>0){
        n.snp.total <- 0
        fam.file <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
        prs.score <- rep(0,nrow(fam.file))
        for(j in 1:22){
          
          #get the number of
          idx <- which(prs.all$CHR==j)
          n.snp.total = n.snp.total+length(idx)
          if(length(idx)>0){
            
            
            filename <- paste0(temp.dir,"prs_chr_",j,".p_value_",k,".profile")
            
            prs.temp <- fread(filename)  
            prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
          }
        }
        write.table(prs.score,file = paste0(out.dir,"/prs_pvalue_",k,".profile"),row.names = F,col.names = F,quote=F)
        
      
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
