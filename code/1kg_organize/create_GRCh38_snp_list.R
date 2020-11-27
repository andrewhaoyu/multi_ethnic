#create grch38 snp list 
setwd("/dcl01/chatterj/data/jzhang2/1000G/GRCh38/EUR")
library(data.table)
#library(bigreadr)
snp.infor.list <- list()
for(i in 1:22){
  print(i)
  temp = fread(paste0("chr",i,".bim"))
  snp.infor.list[[i]] = temp
}
snp.infor.38 <- rbindlist(snp.infor.list)
colnames(snp.infor.38) <- c("CHR","rs_id","unknown","position_GRCh38","REF","ALT")
library(dplyr)
snp.infor.38 = snp.infor.38 %>% 
  select("rs_id","position_GRCh38")
save(snp.infor.38,file = "/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/snp.infor.38.rdata")
