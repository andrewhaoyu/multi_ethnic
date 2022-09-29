data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
summary.eur <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
for(i in 1:5){
  sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[3],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
  bim1 = sum  %>% 
    mutate(V3=0) %>% 
    select(CHR,rsid,V3,BP,A1,A2)
  
  sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[4],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
  bim2 = sum  %>% 
    mutate(V3=0) %>% 
    select(CHR,rsid,V3,BP,A1,A2)
  print(all.equal(bim1,bim2))
  
}
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
prs_cs_ref = read.table("/gpfs/gsfs11/users/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3",header=T)
prs_cs_ref = prs_cs_ref %>% select(SNP)
summary.eur = inner_join(summary.eur,prs_cs_ref,by = c("rs_id"="SNP"))

bim_list = list()
for(i in 1:5){
  bim = read.table(paste0(cur.dir,eth[i],"/all_chr_test.mega.bim"))
  
  bim_list[[i]] = bim  
   
  
}
bim_update = rbindlist(bim_list) %>% distinct(V2, .keep_all = TRUE)
write.table(bim_update, file = paste0(out.dir, "all_eth.bim"),
            row.names = F, col.names = F, quote = F)
