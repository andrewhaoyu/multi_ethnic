#goal: merge ARIC different chrs into one dataset

data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"

i = 1
l = 1
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")

total = 21
temp = 1
filename = rep("c",total)
for(j in 2:22){
  filename[temp] = paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j)
  temp = temp+1
}
out.dir = paste0("/dcl02/leased/chatterj/hzhang1/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
write.table(filename,file = paste0(temp.dir,"merge_list.txt"),row.names = F,col.names = F,quote=F)

system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",1," --merge-list ",temp.dir,"merge_list.txt --make-bed --out ",temp.dir,"all_chr"))

#transform into rds using bigsnpr
library(bigsnpr)
snp_readBed(paste0(temp.dir,"all_chr.bed"))




i = 2
l = 1
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")

total = 21
temp = 1
filename = rep("c",total)
for(j in 2:22){
  filename[temp] = paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j)
  temp = temp+1
}
out.dir = paste0("/dcl02/leased/chatterj/hzhang1/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
write.table(filename,file = paste0(temp.dir,"merge_list.txt"),row.names = F,col.names = F,quote=F)

system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",1," --merge-list ",temp.dir,"merge_list.txt --make-bed --out ",temp.dir,"all_chr"))
snp_readBed(paste0(temp.dir,"all_chr.bed"))
