library(bigreadr)
library(data.table)
library(dplyr)
race_vec = c("EUR","AFR","AMR","EAS")
for(i in 1:5){
  race = race_vec[i]
  sumraw = fread(paste0('/dcs04/nilanjan/data/jjin/tl/',race,'/sumdat/sumdat.txt'))
  if (chr == 6){
    sumraw = sumraw[(sumraw$POS>29e6 & sumraw$POS<33e6),]
  }
  
  valdat_list = list()
  for(chr in 1:22){
    temp_data = fread(paste0('/dcs04/nilanjan/data/jjin/tl/',race,'/geno/chr.qc',chr,'.bim'))
    valdat_list[[chr]] = temp_data
  }
  valdat = rbindlist(valdat_list)
  names = colnames(sumraw)
  sum_dat = inner_join(sumraw,valdat,by=c("SNP_ID"="V2")) %>% 
    select(names)
  save(sum_dat, file = paste0("/dcs04/nilanjan/data/hzhang1/TEL_dat/sum_dat_",race,".rdata"))
  
}
