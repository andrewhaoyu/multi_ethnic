args = commandArgs(trailingOnly = T)
j = as.numeric(args)
#for(j in 1:22){
  print(j)
  system(paste0("zgrep -v \"^##\" /dcl01/chatterj/data/jzhang2/1000G/GRCh37/original_files/ALL.chr",j,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | cut -f1-5 > /dcs04/nilanjan/data/hzhang1/multi_ethnic/result/KG.plink/chr_",j,"_infor"))
  
 # }
  
  library(stringr)
 filename = str_c(paste0("chr_",c(1:22),"_infor "),collapse ="")
  
 system(paste0("cd /dcs04/nilanjan/data/hzhang1/multi_ethnic/result/KG.plink/; head -1 chr_1_infor > chr_head"))
 system(paste0("cd /dcs04/nilanjan/data/hzhang1/multi_ethnic/result/KG.plink/; awk 'FNR>1' ",filename," > tempfile"))

 system(paste0("cd /dcs04/nilanjan/data/hzhang1/multi_ethnic/result/KG.plink/; cat chr_head tempfile > chr_all"))
 
 
 setwd("/data/zhangh24/multi_ethnic/result/LD_simulation_new")
 mega.id <- fread("mega-hm3-rsid.txt",header=F)
 KG.infor <- fread("/data//zhangh24/KG.plink/KG.infor.rsid",header=T)
 
 mega.match = inner_join(mega.id,KG.infor,by=c("V1"="ID"))
 colnames(mega.match)[2] = "CHR"
 save("mega.match",file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/mega.match.rdata"))
 