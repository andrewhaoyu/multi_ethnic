#merge 1kg mega to all_chr
eth = c("EUR", "AFR", "AMR", "EAS", "SAS")
for(i in 1:5){
  data.dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
  soft.dir = "/data/zhangh24/software/"
  file1 = paste0(data.dir,"chr1")
  
  merge_list = paste0(data.dir,"merge_list.txt")
  system(paste0(soft.dir,"plink2 ",
                "--bfile ",file1, " ",
                "--merge-list ",merge_list," ",
                "--make-bed ",
                "--keep-allele-order ",
                "--out ", data.dir, "all_chr"))
}

