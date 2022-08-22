i = 2
m = 1
v = 2
i1 = 1

l = 1
i_rep = 1
library(data.table)
setwd(paste0(out.dir.sum,eth[i],"/prscsx_mega/"))

prs_mega = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi",phi[v],".txt"))

phi = c("1e+00","1e-02","1e-04","1e-06")
setwd(paste0(out.dir.sum,eth[i],"/prscsx/"))
prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi",phi[v],".txt"))


idx = order(abs(prs$V6),decreasing = T)
prs_temp = prs[idx,]



colnames(prs) = c("CHR","SNP","POS","A1_hm3","A2_hm3","effect_hm3")
colnames(prs_mega) = c("CHR","SNP","POS","A1_mega","A2_mega","effect_mega")

prs_mega_temp = prs_mega %>% select(SNP, A1_mega, A2_mega, effect_mega)

idx = order(abs(prs_mega_temp$effect_mega),decreasing = T)
prs_mega_temp = prs_mega_temp[idx,]

idx <- which(prs_mega_temp$SNP%in%prs_temp$SNP==F)
prs_mega_sole = prs_mega_temp[idx,]
prs_com = inner_join(prs_mega_temp,prs_temp,by = "SNP")
head(prs_com)
save(prs_com, file = "/data/zhangh24/multi_ethnic/temp.rdata")
