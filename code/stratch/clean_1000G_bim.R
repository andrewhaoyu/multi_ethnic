#clean the bim file of 1000G

args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])

eth_vec = c("EUR","AFR","AMR","EAS","SAS")
library(data.table)

start_path_vec = c("/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37/",
               "/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh38/",
               "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh37/",
               "/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh38/")
for(l in 1:4){
  start_path = start_path_vec[l]
  
  file_path <- paste0(start_path, 
                      eth_vec[i], "/",
                      "chr", j)
  bim_file_path = paste0(start_path, 
                         eth_vec[i], "/",
                         "chr", j,".bim")
  bim = fread(bim_file_path)
  idx <- which(bim$V2==".")
  # exclude_snp_list = "."
  # write.table(exclude_snp_list, file = "/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37/exclude_snp_list",
  #             row.names  = F, col.names = F, quote = F)
  if(length(idx)>=0){
    res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                        "--threads 2 ",
                        "--bfile ",file_path," ",
                        "--exclude /data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37/exclude_snp_list ",
                        "--make-bed ",
                        "--out ",start_path,eth_vec[i], "/chr",j))
    
  }
  
}
