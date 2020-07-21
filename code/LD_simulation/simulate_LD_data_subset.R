#simulate genotype data by chromosome and by ethnic groups
#EUR,AFR,AME,EAS have 120*k genotype data
#since hapgen2 required at least one disease SNP, we need to know the position information for one SNP
#simulate data for AFR, AMR, EAS
#use /lstractch space to save the output
#then subset the output to save in local space
args <- commandArgs(trailingOnly = T)
#i represent ethnic groups
#j represent chromosome
#k is the simulate replicate
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])
k = as.numeric(args[[3]])
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'))
n <- c(120,120,120,120,120)

tag<- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i],"_chr",j,".tag"),header=F)

system(paste0("/data/zhangh24/software/hapgen2 -m /data/zhangh24/KG.impute2/1000GP_Phase3/genetic_map_chr",j,"_combined_b37.txt -l /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend -h /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap -o /lscratch/",sid,"/test/",eth[i],"_chr",j,"_",k," -n ",n[i]," 1 -dl ",tag[1,1]," 1 1 1 -no_haps_output"))  
  

#hapgen2 simulated the genotype data for all the SNPs
#we only need to genotype data for SNPs in the tag file for each population
#the tag file was created based on MAF> 1% in either population
#tag file only contains the snp position columns, no snpid, chr, allele informatin inside
#tag info file is created for all the snps information

library(data.table)
library(dplyr)
eth <- c("AFR","AMR","EAS","EUR","SAS")
#eth <- c("EUR","AFR","AMR","EAS")


gen <- as.data.frame(fread(paste0("/lscratch/",sid,"/test/",eth[i],"_chr",j,"_",k,".controls.gen"),header=F))
colnames(tag) <- "position"
colnames(gen)[3] <- "position"
tag.gen <- left_join(tag,gen,by="position")
#reorder the column back to the original order
tag.gen <- tag.gen[,c(6:ncol(tag.gen))]
write.table(tag.gen,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,"_",k,".controls.tag.gen"),quote = F,col.names = F,row.names = F)

#we will need to combine tag info file with tag gen file for each ethnic groups



