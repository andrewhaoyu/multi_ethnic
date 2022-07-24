path_to_ref = "/data/zhangh24/test1/ref1/"
path_to_bim = "/data/zhangh24/software/PRScsx/test_data/test"
path_to_sum = "/data/zhangh24/software/PRScsx/test_data/"
out_dir = "/data/zhangh24/software/PRScsx/test_data/"

system(paste0("python /data/zhangh24/software/PRScsx/PRScsx.py", 
              " --ref_dir=",path_to_ref,
              " --bim_prefix=",path_to_bim,
              " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,eth[i],"_sumstats.txt",
              " --n_gwas=200000,100000",
              " --pop=EUR,EAS",
              " --chrom=22",
              " --phi=1e-02", 
              " --out_dir=",out_dir,
              " --out_name=test1"))


path_to_ref = "/data/zhangh24/test1/ref2/"
path_to_bim = "/data/zhangh24/test1/PRScsx/test_data/test"
path_to_sum = "/data/zhangh24/test1/PRScsx/test_data/"
out_dir = "/data/zhangh24/test1/PRScsx/test_data/"

system(paste0("python /data/zhangh24/test1/PRScsx/PRScsx.py -h"))

system(paste0("cd /data/zhangh24/test1/; ",
              "python ./PRScsx/PRScsx.py", 
              " --ref_dir=./ref2",
              " --bim_prefix=./PRScsx/test_data/test",
              " --sst_file=./PRScsx/test_data/EUR_sumstats.txt,./PRScsx/test_data/EAS_sumstats.txt",
              " --n_gwas=200000,100000",
              " --pop=EUR,EAS",
              " --chrom=22",
              " --phi=1e-02", 
              " --out_dir=./",
              " --out_name=test2"))


test1 = fread(paste0(out_dir,"test1_EAS_pst_eff_a1_b0.5_phi1e-02_chr22.txt"))
idx <- which(test1$V2 == "rs174342")
test1[idx,]



test2 = fread(paste0(out_dir,"test2_EAS_pst_eff_a1_b0.5_phi1e-02_chr22.txt"))
idx <- which(test2$V2 == "rs174342")
test2[idx,]
idx <- which(snp_mult$SNP=="rs174342")
snp_mult[idx,]



install.packages("BiocManager")
BiocManager::install("rhdf5")
library(rhdf5)
path_to_ref = "/data/zhangh24/software/PRScsx/1KGLD/"
h5f = H5Fopen(paste0(path_to_ref,"/ldblk_1kg_amr/ldblk_1kg_chr1.hdf5"))
h5f2 = H5Fopen(paste0(path_to_ref,"/ldblk_1kg_afr/ldblk_1kg_chr1.hdf5"))
head(h5f$blk_1)
all_snp_list = list()
for(i in 1:24){
  all_snp_list[[i]] = data.frame(SNP = eval(parse(text = paste0("h5f$blk_",i,"$snplist"))))
  
}

snp_data = data.frame(SNP = eval(parse(text = paste0("h5f$blk_",2,"$snplist"))))

all_snp = rbindlist(all_snp_list)
idx <- which(all_snp$SNP=="rs174342")
all_snp[idx,]





path_to_ref = "/data/zhangh24/test1/ref1/"
h5f = H5Fopen(paste0(path_to_ref,"/ldblk_1kg_eas/ldblk_1kg_chr22.hdf5"))
all_snp_list = list()
for(i in 1:24){
  all_snp_list[[i]] = data.frame(SNP = eval(parse(text = paste0("h5f$blk_",1,"$snplist"))))
  
}

all_snp = rbindlist(all_snp_list)
idx <- which(all_snp$SNP=="rs174342")
all_snp[idx,]

