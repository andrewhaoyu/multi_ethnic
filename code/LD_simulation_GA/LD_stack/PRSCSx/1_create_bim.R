data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
eth <- c("EUR","AFR","AMR","EAS","SAS")

bim_list = list()
for(i in 1:5){
  bim = read.table(paste0(cur.dir,eth[i],"/all_chr_test.mega.bim"))
  
  bim_list[[i]] = bim  
   
  
}
bim_update = rbindlist(bim_list) %>% distinct(V2, .keep_all = TRUE)
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")

bim_update = as.data.frame(fread(paste0(out.dir, "all_eth.bim")))

snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)

bim_update = left_join(bim_update,snp.infor,by = c("V2"="SNP"))

bim_update = bim_update %>% 
  select(V1, rs_id, V3, V4, V5, V6) %>% na.omit()

write.table(bim_update, file = paste0(out.dir, "all_eth.bim"),
            row.names = F, col.names = F, quote = F)
