path_to_ref = "/data/zhangh24/test1/ref1/"
path_to_bim = "/data/zhangh24/software/PRScsx/test_data/test"
path_to_sum = "/data/zhangh24/software/PRScsx/test_data/"
out_dir = "/data/zhangh24/software/PRScsx/test_data/"

system(paste0("python /data/zhangh24/software/PRScsx/PRScsx.py", 
              " --ref_dir=",path_to_ref,
              " --bim_prefix=",path_to_bim,
              " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,"EAS_sumstats.txt",
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

#system(paste0("python /data/zhangh24/test1/PRScsx/PRScsx.py -h"))

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
path_to_ref = "/data/zhangh24/test1/ref2/"
snp_mult = fread(paste0(path_to_ref,"snpinfo_mult_1kg_hm3"))
idx <- which(snp_mult$SNP=="rs1001276")
h5f = H5Fopen(paste0(path_to_ref,"/ldblk_1kg_eur/ldblk_1kg_chr1.hdf5"))

 data.frame(SNP = h5f$blk_1$snplist)


for(i in 1:133){
  snp_list = eval(parse(text = paste0("h5f$blk_",52)))
  data  = data.frame(SNP = eval(parse(text = paste0("h5f$blk_",52,"$snplist"))))
  #data2  = data.frame(SNP = eval(parse(text = paste0("h5f$blk_",6,"$snplist"))))
  idx <- which(data$SNP%in%"rs1001276")
  LD = data.frame(SNP = eval(parse(text = paste0("h5f$blk_",52,"[[1]]"))))
  #idx2 = which(data2$SNP%in%"rs1001276")
  
  if("rs1001276"%in%data$SNP==T){
    print(i)
  }
}

 data[846:849,]

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

