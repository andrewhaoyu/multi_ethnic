args = commandArgs(trailingOnly = T)
j = as.numeric(args[[1]])


#remove all the multi-allelic snps
library(data.table)
library(dplyr)
system(paste0("/data/zhangh24/software/plink2 --vcf /gpfs/gsfs10/users/BC_risk_prediction/ghana/GSA/VCFs_to_share/subtracted_chr",j,".dose.vcf.gz --recode --make-bed --out /data/zhangh24/multi_ethnic/data/GBHS_plink/chr",j," --maf 0.001 --double-id"))
bim <- fread(paste0("/gpfs/gsfs11/users/zhangh24/multi_ethnic/data/GBHS_plink/chr",j,".bim"),header=F)
id = bim %>% mutate(id = paste0(V2,":",V6,":",V5)) %>% 
  select(id)
bim$V2 = id
write.table(bim,paste0("/gpfs/gsfs11/users/zhangh24/multi_ethnic/data/GBHS_plink/chr",j,".bim"),row.names = F,col.names = F,quote = F)
# dup.id <- bim$V2[duplicated(bim$V2)]
# write.table(dup.id, file = paste0("/data/zhangh24/multi_ethnic/data/GHBS_plink/dup.id.chr.",j),row.names = F,col.names = F,quote=F)
# system(paste0("/data/zhangh24/software/plink2 --vcf /gpfs/gsfs10/users/BC_risk_prediction/ghana/GSA/VCFs_to_share/subtracted_chr",j,".dose.vcf.gz --recode --make-bed --out /data/zhangh24/multi_ethnic/data/GHBS_plink/chr",j," --maf 0.005 --double-id --biallelic-only --exclude /data/zhangh24/multi_ethnic/data/GHBS_plink/dup.id.chr.",j))
# 


# merge all the dup id
# dup.id.list = list()
# for(j in 1:22){
#  dup.id <- fread(paste0("/data/zhangh24/multi_ethnic/data/GHBS_plink/dup.id.chr.",j),header=F)
#  dup.id.list[[j]] = dup.id
# }
# dup.id = rbindlist(dup.id.list)
# 
# save(dup.id,file = paste0("/data/zhangh24/multi_ethnic/data/GHBS_plink/all.dup.id.rdata"))








# for(j in 1:22){
# 
# }

#dup.id <- fread(paste0("/data/zhangh24/multi_ethnic/data/GHBS_plink/dup.id.chr.",j),header=F)

# for(j in 1:22){
#   system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/multi_ethnic/data/GHBS_plink/chr",j," --exclude /data/zhangh24/multi_ethnic/data/GHBS_plink/all_chr-merge.missnp --make-bed --out /data/zhangh24/multi_ethnic/data/GHBS_plink/chr_clean",j))
# }



# dup.id <- bim$V2[duplicated(bim$V2)]
# 
# idx <- which(bim$V2=="1:15274")
# /data/zhangh24/software/bcftools norm -m-any --check-ref /gpfs/gsfs10/users/BC_risk_prediction/ghana/GSA/VCFs_to_share/subtracted_chr22.dose.vcf.gz | /data/zhangh24/software/bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' > /data/zhangh24/multi_ethnic/data/GHBS_plink/chr22.bcf
#   
# 
# 
# 
# 
# annotate -- /gpfs/gsfs10/users/BC_risk_prediction/ghana/GSA/VCFs_to_share/subtracted_chr22.dose.vcf.gz