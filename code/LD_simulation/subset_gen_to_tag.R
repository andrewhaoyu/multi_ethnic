#hapgen2 simulated the genotype data for all the SNPs
#we only need to genotype data for SNPs in the tag file for each population
#the tag file was created based on MAF> 1% in either population
#tag file only contains the snp position columns, no snpid, chr, allele informatin inside
#tag info file is created for all the snps information

library(data.table)

#eth <- c("AFR","AMR","EAS","EUR","SAS")
eth <- c("EUR","AFR","AMR","EAS")
tag_all <- list()
args = commandArgs(trailingOnly = T)
#i1 represent ethnic groups
#i2 represent chromosome
#i3 is the simulate replicate
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
i3 = as.numeric(args[[3]])
#tag file contains the SNPs information with MAF>=0.01 for specific ethnic

tag <- read.table(paste0("/data/zhangh24/KG.impute2/tag/",eth[i1],"_chr",i2,".tag"),header=F)
  
library(dplyr)
    
      gen <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",i2,"_",i3,".controls.gen"),header=F))
      colnames(tag) <- "position"
      colnames(gen)[3] <- "position"
      tag.gen <- left_join(tag,gen,by="position")
      #reorder the column back to the original order
      tag.gen <- tag.gen[,c(6:ncol(tag.gen))]
      write.table(tag.gen,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i1],"/chr",i2,"_",i3,".controls.tag.gen"),quote = F,col.names = F,row.names = F)
    


