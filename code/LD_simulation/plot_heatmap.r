#prepare the sample data for different ethnic groups for ploting
#cd /data/zhangh24/multi_ethnic/result/LD_simulation/EUR
#head -100  chr22.controls.tags.gen > chr22.controls.tags_100.gen
#similar process for all the ethnic groups

#use qctools to tranform the 1000 genome haplotype data into the .gen format
#only do this only chr 1
#need to have the legend file for every ethnic groups
#cp 1000GP_Phase3_chr1.legend ../EUR/1000GP_Phase3_chr1.legend 
#need to change the hap and legend name to be the same
#cd /data/zhangh24/KG.impute2/EUR
#mv 1000GP_Phase3_chr1.legend chr1.legend
#same procedure for all other ethnic groups

code <- rep("c",4)
eth <- c("EUR","AFR","AMR","EAS")
for(k in 1:4){
  code[k] <- paste0("/data/zhangh24/software/qctool/qctool -g /data/zhangh24/KG.impute2/",eth[k],"/chr1.hap -filetype impute_haplotypes -og /data/zhangh24/KG.impute2/",eth[k],"/1000GP_Phase3_MTHFR.gen -incl-positions /data/zhangh24/KG.impute2/tag/MTHFR.tag")
  
}

write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/transform_1kg_hap_to_gen.sh"),row.names = F,
            col.names = F,
            quote=F)


#download 1kg data and simulated data to local /Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation/
#load data for plot
library(corrplot)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
setwd("/Users/zhangh24/GoogleDrive")
for(i in 1:4){
  data <- as.data.frame(fread(paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/1000GP_Phase3_MTHFR.gen")))
  
  #impute.form every 3 rows represent one subject
  #first row aa, second row Aa, third row AA
  impute.form <- t(data[,-c(1:6)])
  n.sub <- nrow(impute.form)/3
  n.snp <- ncol(impute.form)
  genotype <- matrix(0,n.sub,n.snp)
  for(k in 1:n.sub){
    genotype[k,] <- impute.form[(3*k-1),]+2*impute.form[3*k,]
      
  }
  MAF <- colSums(genotype)/(2*n.sub)
  idx <- which(MAF>=0.01&
                 MAF<=0.99 )
  genotype <- genotype[,idx]
  LD <- cor(genotype)^2
  png(file = paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/1kg_heatmap.png"),width=6,height=6,units="in",res=300)
  corrplot(LD,method="color",type="upper",diag = T,tl.cex = 1,tl.col = "white")
  dev.off()
  data2 <- as.data.frame(fread(paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/chr1_MTHFR.gen")))
  impute.form <- t(data2[,-c(1:5)])
  n.sub <- nrow(impute.form)/3
  n.snp <- ncol(impute.form)
  genotype <- matrix(0,n.sub,n.snp)
  for(k in 1:n.sub){
    genotype[k,] <- impute.form[(3*k-1),]+2*impute.form[3*k,]
    
  }
  
  for(k in 1:n.sub){
    genotype[k,] <- impute.form[(3*k-1),]+2*impute.form[3*k,]
    
  }
  MAF <- colSums(genotype)/(2*n.sub)
  idx <- which(MAF>=0.01&
                 MAF<=0.99 )
  genotype <- genotype[,idx]
  LD2 <- cor(genotype)^2
  png(file = paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/simulated_heatmap.png"),width=6,height=6,units="in",res=300)
  corrplot(LD2,method="color",type="upper",diag = T,tl.cex = 1,tl.col = "white")
  dev.off()

  }
