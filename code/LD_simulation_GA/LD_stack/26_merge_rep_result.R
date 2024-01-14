#merge each rep result for simulation
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
load(paste0(out.dir,"LD.clump.result.CT.rep.rdata"))
result.table = result.table[,c(1,2,3,5,4,6,7)]
result_list = list()
names = c("eth_vec",
          "r2_vad",
          "l_vec",
          "m_vec",
          "ga_vec",
          "rep_vec",
          "method_vec")
#ct
ct.result = result.table
colnames(ct.result) = names
result_list[[1]] = result.table
#eur ct
load(paste0(out.dir,
       "eur.snp.reult.rep.rdata"))
colnames(eursnp.result) = names
result_list[[2]] = eursnp.result
#weigted prs
load(paste0(out.dir,"weightedprs.CT.rep.rdata"))
weighted.ct = result.table
colnames(weighted.ct) = names
result_list[[3]] = weighted.ct
#prs-csx
load(paste0(out.dir,
       "prscsx_five.result_rep_ga.rdata"))
colnames(prscsx.result) = names
result_list[[4]] = prscsx.result
#ctsleb
load(paste0(out.dir,"LD.clump.result.allethtest.rep.rdata"))
colnames(alleth.EB.result.rep) = names
result_list[[5]] = alleth.EB.result.rep
result_rep = rbindlist(result_list)
save(result_rep,file= paste0(out.dir,"simu_result_rep.rdata"))


result_rep %>%  filter(l_vec==3&m_vec==3&ga_vec==1)
library(dplyr)
result_sub = result_rep %>% filter(ga_vec==5)
table(result_sub$method_vec)

load("./result/LD_simulation_GA/LD_stack/simu_result_rep.rdata")


