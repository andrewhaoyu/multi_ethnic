#goal: create .ld file using plink
#goal: create block_chr, block_size and snp_list for prs_csx
#prscsx requires three files to create hdf file
#block_chr is the chromosome number for each block
#block_size is the number of SNPs in LD block
#snp_list is the SNPs in each LD block
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
library(bigsnpr)
library(data.table)
library(dplyr)
library(rlang)
#load LD blocks table from Berisa, T.  Bioinformatics, 2016.
eth = c("EUR","AFR","AMR","EAS","SAS")
pos_table = fread(paste0("/data/zhangh24/MR_MA/data/ld_block"))


#data dir for 1KG reference with MEGA chip
data_dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
#out_dir saves the ld output from plink
out_dir = paste0("/data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"/block_ld/")
#snp_list saves the block_chr, block_size and snp_list for write_hdf
snp_list_dir = paste0(out_dir,"snplist_ldblk/")
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp_dir = paste0('/lscratch/',sid,'/test/')
#generate SNP list for each block


pos_table_sub = pos_table %>% 
  as.data.frame() %>% 
  mutate(CHR_clean = as.numeric(gsub("chr","",x = chr)))

#get the start index and end index for each chromosome
start = rep(0, 22)
end = rep(0,22)
for(j in 1:22){
  idx <- which(pos_table_sub$CHR_clean==j)
  start[j] = idx[1]
  end[j] = idx[length(idx)]
}


load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.match37_38.rdata")


snp.infor.match.select = snp.infor.match %>% 
  select(id,rs_id,AFR,AMR,EAS,EUR,SAS,position,CHR)

bim_file = fread(paste0(data_dir,"all_chr.bim"))



bim_file = left_join(bim_file,snp.infor.match.select,by=c("V2"="rs_id"))

#block_chr is the chromosome number for each block
block_chr = pos_table_sub %>% 
  mutate(CHR_clean = as.numeric(gsub("chr","",x = chr))) %>% 
  select(CHR_clean)
write.table(block_chr,file = paste0(snp_list_dir,"blk_chr"),row.names = F,col.names = F,quote = F)

#block_size is the number of SNPs in LD block



block_size = data.frame(size = rep(0,nrow(pos_table_sub)))


soft_dir = "/data/zhangh24/software/"



for(k in 1:nrow(pos_table_sub)){
  #if it's first block, increase the block size to include SNPs before start
  if(k %in% start ==1){
    start =1
    end = pos_table_sub[k,3]
  }else if(k %in% end ==1){
    #if it's last block, increase end to include any SNPs
    start = pos_table_sub[k,2]
    end = Inf
  }else{
    start = pos_table_sub[k,2]
    end = pos_table_sub[k,3]
  }
  
  j = block_chr[k,1]
  print(j)
  #subset to common SNPs in the population
  #subset to bi-allelic snps
  snp_list = bim_file %>% 
    filter(CHR == j) %>% 
    filter(V4>=start&V4<=end) %>% 
    filter(get(eth[i])>=0.01&
             get(eth[i])<=0.99) %>% 
    filter(V5 %in% c("A", "T", "C", "G")&
             V6 %in% c("A", "T", "C", "G")) %>%
    select(V2) %>% 
    rename(SNP=V2)
  
  if(nrow(snp_list)!=0){
    write.table(snp_list,file= paste0(temp_dir,"extract_snp_list"),row.names = F,col.names = T,quote=F)  
    #snp_list is the SNPs in each LD block
    
    snp_list_vec = snp_list[l,1]
    for(l in 2:nrow(snp_list)){
      
        snp_list_vec = paste0(snp_list_vec," ",snp_list[l,1])  
      
    }
    write.table(snp_list_vec,file= paste0(snp_list_dir,"snplist_blk",k),row.names = F,col.names = F,quote=F)
    block_size[k,1] = nrow(snp_list)
    #write.table(snp_list,file= paste0(out_dir,"snp_list_100"),row.names = F,col.names = T,quote=F)
    #calculate LD for selected block
    # system(paste0(soft_dir,"plink2 ",
    #               "--bfile ",data_dir,"chr",j," ",
    #               "--keep-allele-order ",
    #               "--extract ",temp_dir,"extract_snp_list ",
    #               "--r square ",
    #               "--out ", temp_dir,"ldblock_",k))
    #system(paste0("more /lscratch/42052892/test/ldblock_100.ld"))
    # system(paste0("mv ", temp_dir,"ldblock_",k,".ld ",
    #               out_dir))
    # 
  }else{
    block_size[k,1] = nrow(snp_list)
  }
  
  
  
  

}
write.table(block_size, file= paste0(snp_list_dir,"blk_size"),row.names = F,col.names = F,quote=F)

