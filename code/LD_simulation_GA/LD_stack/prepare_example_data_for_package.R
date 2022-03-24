#goal: prepare example data for R package
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
load(paste0(out.dir.sum,"time_mem/eur_sumdata.rdata"))
library(dplyr)
library(data.table)
summary.eur.select = summary.eur %>% 
  mutate(SE = BETA/STAT) %>% 
  rename(P=peur) %>% 
  select(CHR,SNP,BP,A1,BETA,SE,P,rs_id)
write.table(summary.eur.select,
            paste0(out.dir.sum,"time_mem/EUR_sumdata.txt"),
            quote=F,
            row.names = F,
            col.names = T)
load(paste0(out.dir.sum,"time_mem/tar_sumdata.rdata"))
summary.tar.select = summary.tar %>% 
  mutate(SE = BETA/STAT) %>% 
  #rename(P=peur) %>% 
  select(CHR,SNP,BP,A1,BETA,SE,P,rs_id)
write.table(summary.tar.select,
            paste0(out.dir.sum,"time_mem/AFR_sumdata.txt"),
            quote=F,
            row.names = F,
            col.names = T)

i = 2
l = 3
m = 1
i_rep = 1
i1 = 1
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
colnames(mega.list) = "rs_id"
eth <- c("EUR","AFR","AMR","EAS","SAS")
lower_eth = tolower(eth[i])
for(i in 3:5){
  summary <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
  summary.match = left_join(summary,snp.infor,by="SNP")
  summary.match.update = inner_join(summary.match,mega.list,by="rs_id")
  summary.select = summary.match.update %>% 
    filter(CHR==22) %>% 
    mutate(SE = BETA/STAT) %>% 
    #rename(P=peur) %>% 
    select(CHR,SNP,BP,A1,BETA,SE,P,rs_id)
  write.table(summary.select,
              paste0(out.dir.sum,"time_mem/",eth[i],"_sumdata.txt"),
              quote=F,
              row.names = F,
              col.names = T)
}


out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
load(paste0(out.dir,"y_test.rdata"))
load(paste0(out.dir,"y_vad.rdata"))
write.table(y_test,
            paste0(out.dir.sum,"time_mem/","y_tuning.txt"),
            quote=F,
            row.names = F,
            col.names = F)

write.table(y_vad,
            paste0(out.dir.sum,"time_mem/","y_validation.txt"),
            quote=F,
            row.names = F,
            col.names = F)