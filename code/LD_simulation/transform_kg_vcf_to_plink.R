#LD clumping require process vcf (hap,lengend) to plink
#use qctool version 2 to process it
#merge all the 22 chr vcf files to one file for each ethnic group
eth <- c("EUR","AFR","AMR","EAS","SAS")
code <- rep("c",1000)
temp <- 1
#merge hap files together by ethnic groups
for(i in 1:5){
  temp.code <- paste0("cat ")
  for(j in 1:22){
   temp.code <- paste0(temp.code," /data/zhangh24/KG.impute2/",eth[i],"/chr",j,".hap") 
  }
  temp.code <- paste0(temp.code, " > /data/zhangh24/KG.impute2/",eth[i],"/chr_all.hap")
  code[temp] <- temp.code
  temp <- temp+1
}
#merge legend file together
#don't need to be ethnic specific
#since the hap files contain all the SNPs
#if a specific SNP doesn't exist in a particular ethnic group
#then the row of the hap file will all be 0

  temp.code <- paste0("head -1 /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr1.legend > /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr_all.legend; awk 'FNR>1' ")
  for(j in 1:22){
    temp.code <- paste0(temp.code," /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",j,".legend") 
  }
  temp.code <- paste0(temp.code, " >> /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr_all.legend")
  code[temp] <- temp.code
  temp <- temp+1

code <- code[1:(temp-1)]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/merge_hap_legend_together.sh",row.names = F,col.names = F,quote=F)

#move the 1000GP_Phase3_chr_all.legend to different ethnic folder to prepare for gtool to transform vcf to plink
#cp /data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr_all.legend /data/zhangh24/KG.impute2/EUR/chr_all.legend



#use gtool to transform vcf (hap,legend) to plink
code <- rep("c",1000)
temp <- 1
for(i in 5:5){
    code[temp] <- paste0("/data/zhangh24/software/qctool/qctool -g /data/zhangh24/KG.impute2/",eth[i],"/chr_all.hap -filetype impute_haplotypes -og /data/zhangh24/KG.plink/",eth[i],"/chr_all  -ofiletype binary_ped")
    temp <- temp+1
  
}
code <- code[1:(temp-1)]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/transform_1kg_hap_to_plink.sh",row.names = F,col.names = F,quote=F)


#add the chromosome information to the bim file
library(data.table)
library(dplyr)
bim <- as.data.frame(fread("/data/zhangh24/KG.plink/EUR/chr_all.bim",header=F))
eth <- c("EUR","AFR","AMR","EAS")
CHR <- rep(0,nrow(bim))
total <- 0
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
 temp <- nrow(leg)
 CHR[total+(1:temp)] <- i
 total <- total+temp
    
  }
bim[,1] <- CHR  
write.table(bim,file = "/data/zhangh24/KG.plink/EUR/chr_all.bim",row.names = F,col.names = F,quote=F)
