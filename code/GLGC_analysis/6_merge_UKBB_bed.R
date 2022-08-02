#goal: merge UKBB data into a single file
#create merge-list file
eth = c("EUR", "AFR", "AMR", "EAS", "SAS")


for(i in 1:5){
  data.dir = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/",eth[i],"/")
  file = data.frame(names = rep("c",21))
  for(j in 2:22){
    file[j-1,1] = paste0(data.dir,"chr",j)
  }
  write.table(file, file = paste0(data.dir, "merge_list.txt"), row.names = F, col.names = F, quote = F)
}

for(i in 1:5){
  data.dir = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/",eth[i],"/")
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
