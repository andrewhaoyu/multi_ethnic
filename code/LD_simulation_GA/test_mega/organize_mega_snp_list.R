library(data.table)
library(dplyr)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
#load mega snp information
mega_snp <- as.data.frame(fread(paste0(cur.dir,"./snpBatch_ILLUMINA_1062317")))
colnames(mega_snp)[5] <- "rs_id"
mega_snp_select = mega_snp %>% 
  select("rs_id","chr","chr_pos") %>% 
  rename(pos = chr_pos) %>% 
  mutate(chr_pos = paste0(chr,":",pos)) %>% 
select(rs_id,chr_pos)
write.table(mega_snp_select,file =paste0(cur.dir,"./mega_snp_chr_pos.txt"),row.names=F,col.names=F,quote=F)
range(as.numeric(mega_snp_select$chr_pos),na.rm = T)

idx <- which(mega_snp_select$chr_pos=="N.D.")

load(paste0(cur.dir,"snp.infor.rdata"))
head(snp.infor)
library(tidyr)
snp.infor_update = snp.infor %>% 
  separate(col=id,into=c("rs_id","position_new","allele1","allele2"),sep=":",remove=F)
  

idx <- which(mega_snp_select$rs_id%in%snp.infor_update$rs_id==F)

jdx <- which(snp.infor_update$position==1335598)

snp_infor_filter
idx <- which(mega_snp$chr_pos==51239678)

mega_hap <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header=F),)
idx <- which(mega_hap$V1=="rs74045277")
inter.id <- intersect(snp.infor_update$rs_id,as.character(mega_hap$V1))
length(inter.id)
