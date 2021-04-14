args = commandArgs(trailingOnly = T)
j = as.numeric(args[[1]])
system(paste0("/data/zhangh24/software/plink2 --vcf /gpfs/gsfs10/users/BC_risk_prediction/ghana/GSA/VCFs_to_share/subtracted_chr",j,".dose.vcf.gz --recode --make-bed --out /data/zhangh24/multi_ethnic/data/GHBS_plink/chr",j," --maf 0.005 --double-id --biallelic-only"))


