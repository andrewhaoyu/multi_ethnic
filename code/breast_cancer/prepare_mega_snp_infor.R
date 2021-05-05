args = commandArgs(trailingOnly = T)
j = as.numeric(args)
#for(j in 1:22){
  print(j)
  system(paste0("zgrep -v \"^##\" /dcl01/chatterj/data/jzhang2/1000G/GRCh37/original_files/ALL.chr",j,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | cut -f1-5 > /dcs04/nilanjan/data/hzhang1/multi_ethnic/result/KG.plink/chr_",j,"_infor"))
  
 # }