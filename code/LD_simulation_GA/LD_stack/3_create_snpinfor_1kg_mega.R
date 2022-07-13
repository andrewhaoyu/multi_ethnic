#goal: snpinfor and snpinfo_mult_1kg_mega for prs_csx
#load LD blocks table from Berisa, T.  Bioinformatics, 2016.
library(data.table)
library(dplyr)
eth = c("EUR","AFR","AMR","EAS","SAS")
pos_table = fread(paste0("/data/zhangh24/MR_MA/data/ld_block"))

#data dir for 1KG reference with MEGA chip
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.match37_38.rdata")


snp.infor.match.select = snp.infor.match %>% 
  select(id,rs_id,AFR,AMR,EAS,EUR,SAS,position,CHR)

#generate SNP infor file for each ethnic group
#bim file for all ancestries are the same
i = 1
data_dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
bim_file = fread(paste0(data_dir,"all_chr.bim"))
bim_file = left_join(bim_file,snp.infor.match.select,by=c("V2"="rs_id"))
snp_info_list = list()
for(i in 1:5){
  
  
  data_dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
  #out_dir saves the ld output from plink
  out_dir = paste0("/data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"/")
  
  snp_info = bim_file %>% 
    filter(get(eth[i])>=0.01&
             get(eth[i])<=0.99) %>%
    select(V1, V2, V4, V5, V6, eth[i]) %>% 
    rename(CHR = V1,
           SNP = V2,
           BP = V4,
           A1  = V5,
           A2 = V6,
           MAF = eth[i])
  snp_info_list[[i]]  = snp_info
  write.table(snp_info, file = paste0(out_dir,"snpinfo_1kg_mega"), row.names = F, col.names = T, quote= F, sep = "\t")
  
}

snp_info_combine = rbindlist(snp_info_list)

unique_snp = data.frame(SNP = unique(snp_info_combine$SNP))


#create snpinfo_mult_1kg_mega
#bim file for all ancestries are the same
#for SNPs not exist in a particular ancestry, the values for the SNP will be 0
#just use EUR to create the reference
i = 1
data_dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
bim_file = fread(paste0(data_dir,"all_chr.bim"))
head(bim_file)
#subset to SNPs with MAF > 0.01 in at least one ancestries
bim_file = left_join(unique_snp,bim_file,by = c("SNP"="V2"))
#order file by CHR and pos
bim_file_update = bim_file[order(bim_file$V1,bim_file$V4),]


bim_file = left_join(bim_file_update,snp.infor.match.select,by=c("SNP"="rs_id"))

FLIP_data = data.frame(matrix(0,nrow(bim_file),length(eth)))
#loop through all ancestries
for(i in 1:5){
bim_file_merge = left_join(bim_file,snp_info_list[[i]],by = "SNP")

FLIP_temp = bim_file_merge %>% 
  mutate(FLIP = case_when(is.na(A1)==1 ~ 0,
                          A1 == V5 ~ 1,
                          A1 == V6 ~ -1)) %>% 
  select(FLIP)
FLIP_data[,i] = FLIP_temp
}
colnames(FLIP_data) = paste0("FLIP_",eth)

bim_file_com = cbind(bim_file, FLIP_data)


snpinfo_mult_1kg = bim_file_com %>% 
  rename(BP = V4,
         A1  = V5,
         A2 = V6,
         FRQ_AFR = AFR,
         FRQ_AMR = AMR,
         FRQ_EAS = EAS,
         FRQ_EUR = EUR,
         FRQ_SAS = SAS) %>% 
  select(CHR, SNP, BP, A1, A2, FRQ_AFR, FRQ_AMR, FRQ_EAS, FRQ_EUR, FRQ_SAS,
         FLIP_AFR, FLIP_AMR, FLIP_EAS, FLIP_EUR, FLIP_SAS)

write.table(snpinfo_mult_1kg, 
            file = paste0("/data/zhangh24/software/PRScsx/1KGLD_MEGA/snpinfo_mult_1kg_mega"), 
            row.names = F, col.names = T, quote= F, sep = "\t")

