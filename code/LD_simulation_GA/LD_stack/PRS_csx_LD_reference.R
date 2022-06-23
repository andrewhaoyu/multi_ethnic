library(bigsnpr)
library(data.table)
library(dplyr)
#load LD blocks table from Berisa, T.  Bioinformatics, 2016.
eth = c("EUR","AFR","AMR","EAS","SAS")
pos_table = fread(paste0("/data/zhangh24/MR_MA/data/ld_block"))
i = 1
j = 1
#data dir for 1KG reference with MEGA chip
data_dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")

out_dir = paste0("/data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"/")
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp_dir = paste0('/lscratch/',sid,'/test/')
#generate SNP list for each block
pos_table_sub = pos_table %>% 
  filter(chr==paste0("chr",j)) %>% 
  as.data.frame()

bim_file = fread(paste0(data_dir,"chr",j,".bim"))
soft_dir = "/data/zhangh24/software/"
for(k in 1:nrow(pos_table_sub)){
  #if it's first block, increase the block size to include SNPs before start
  if(k==1){
    start =1
    end = pos_table_sub[k,3]
  }else if(k==nrow(pos_table_sub)){
    #if it's last block, increase end to include any SNPs
    start = pos_table_sub[k,2]
    end = Inf
  }else{
    start = pos_table_sub[k,2]
    end = pos_table_sub[k,3]
  }
  snp_list = bim_file %>% 
    filter(V4>=start&V4<=end) %>% 
    select(V2) %>% 
    rename(SNP=V2)
  write.table(snp_list,file= paste0(temp_dir,"extract_snp_list"),row.names = F,col.names = T,quote=F)
  #write.table(snp_list,file= paste0(out_dir,"snp_list_100"),row.names = F,col.names = T,quote=F)
  #calculate LD for selected block
  system(paste0(soft_dir,"plink2 ",
                "--bfile ",data_dir,"chr",j," ",
                "--keep-allele-order ",
                "--extract ",temp_dir,"extract_snp_list ",
                "--r square ",
                "--out ", temp_dir,"ldblock_",k))
  system(paste0("more /lscratch/42052892/test/ldblock_100.ld"))
  system(paste0("mv ", temp_dir,"ldblock_",k,".ld ",
                out_dir))
 

}

