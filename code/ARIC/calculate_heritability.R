args = commandArgs(trailingOnly = T)

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
grm.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[1],"/",eth[i],"/")
data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
colnames(y)[2] <- "ID"
fam <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
colnames(fam)[2] <- "ID"
pheno = left_join(fam,y,by="ID") %>% 
  select(V1,ID,y)
write.table(pheno,file =paste0(temp.dir,"pheno.txt"),row.names = F,col.names = F,quote=F)
res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/gcta/gcta64 --grm ",grm.dir,"merged_chr_qc --pheno ",temp.dir,"pheno.txt --reml --out ",out.dir,"heritability"))  
