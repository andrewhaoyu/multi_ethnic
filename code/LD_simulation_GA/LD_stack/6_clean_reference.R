#clean the hd5 reference with nan

eth = c("EUR","AFR","AMR","EAS","SAS")
library(rhdf5)
path_to_ref = "/data/zhangh24/software/PRScsx/1KGLD_MEGA/"
snp_mult = fread(paste0(path_to_ref,"snpinfo_mult_1kg_hm3"))
snp_id_list = list()
id_index = list()
temp = 1
for(i in 1:5){
  pos_table = fread(paste0("/data/zhangh24/MR_MA/data/ld_block/ld_block_",eth[i]))
  colnames(pos_table) = c("chr", "start", "end")
  
  pos_table_sub = pos_table %>% 
    as.data.frame() %>% 
    mutate(CHR_clean = as.numeric(gsub("chr","",x = chr))) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end))
  
  for(j in 1:22){
    h5f = H5Fopen(paste0(path_to_ref,"/ldblk_1kg_",tolower(eth[i]),"/ldblk_1kg_chr",j,".hdf5"))
    snp_set = pos_table_sub %>% filter(CHR_clean==j)
    n_snp_set = nrow(snp_set)
    for(k in 1:n_snp_set){
    print(k)
      
      LD = eval(parse(text = paste0("h5f$blk_",k,"[[1]]")))
      #if this block has snp, then test whether all SNPs LD are clear
      if(nrow(LD)>0){
        row_sum = rowSums(is.nan(LD))  
        idx <- which(row_sum==nrow(LD))
        if(length(idx)>0){
          print(c(i, j, k))
          id_index[[temp]] = data.frame(i_index = i,j_index = j,k_index = k)
          snp_list = eval(parse(text = paste0("h5f$blk_",k,"$snplist")))
          snp_id_list[[temp]] = data.frame(SNP = snp_list[idx])
          temp = temp + 1
        }
      }
      
      
    }  
  }
  
}

snp_id_filter = rbindlist(snp_id_list)
save(snp_id_filter,file = "/data/zhangh24/software/PRScsx/snp_id_filter.rdata")
