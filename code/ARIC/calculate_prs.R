#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is trait
#j is the number of chromsome
args = commandArgs(trailingOnly = T)

j = as.numeric(args[[1]])
#l = as.numeric(args[[3]])
#m = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
#i_rep = 1

#j = as.numeric(args[[3]])
library(dplyr)
library(data.table)

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

#generate q_range file is a one time job
# for(i in 1:2){
#   for(l in 1:3){
#     temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
#     q_range = data.frame(
#       paste0("p_value_",c(1:n_pthres)),
#       rep(0,n_pthres),
#       pthres,stringsAsFactors = F)
#     write.table(q_range,file = paste0(temp.dir,"2Dq_range_file"),row.names = F,col.names = F,quote=F)
#     
#   }
#   
# }

for(i in 1:2){
  for(l in 1:3){
    setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    LD <- as.data.frame(fread(paste0(out.dir,"/LD_clump_chr_",j,".clumped")))
    clump.snp <- LD[,3,drop=F] 
    
    sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    sum.data = sum.data %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL) 
    #load bim file to match all the SNPs with different strand
    #in summary data and ARIC data
    bim <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j,".bim")))
    colnames(bim)[2] <- "SNP"
    sum.com.match <- left_join(bim,sum.data,by="SNP")
    #A1 is in summary data
    #V5 and V6 is from genotype data
    #if A1 is not either V5 or V6, flip strand
    sum.com.match = sum.com.match %>% 
      mutate(A1 = case_when(
        (A1!=V5&A1!=V6)&A1=="A" ~"T",
        (A1!=V5&A1!=V6)&A1=="C" ~"G",
        (A1!=V5&A1!=V6)&A1=="T" ~"A",
        (A1!=V5&A1!=V6)&A1=="G" ~"C",
        TRUE ~ A1
      ))
    
    
    # write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
    #             ,col.names = T,row.names = F,quote=F) 
    #prepare association file for plink
    prs.all <- left_join(clump.snp,sum.com.match,by="SNP")
    n_pthres <- length(pthres)
    #q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres),stringsAsFactors = F)
    prs.file <- prs.all %>% filter(CHR==j) 
    prs.file = prs.file[,c("SNP","A1","BETA")]
    #setwd(temp.dir)
    write.table(prs.file,file = paste0(temp.dir,"prs_coeff_chr_",j),col.names = T,row.names = F,quote=F)
    p.value.file <- prs.all %>% filter(CHR==j) 
    p.value.file = p.value.file[,c("SNP","P")]
      
    # com.prs = left_join(prs.file,p.value.file,by="SNP")
    # com.prs.filter = com.prs %>%
    # filter(P<=pthres[9])

     write.table(p.value.file,file = paste0(temp.dir,"p_value_chr_",j),col.names = T,row.names = F,quote=F)
      
      if(nrow(prs.file)>0){
        res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_chr_",j," header --threads 2 --score ",temp.dir,"prs_coeff_chr_",j," header no-sum no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"prs_chr_",j))
        #res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_chr_",j," header --threads 2 --score ",temp.dir,"prs_coeff_chr_",j," header no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"prs_chr_",j))
        print("step2 finished")
          #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
        if(res==2){
          stop()
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
