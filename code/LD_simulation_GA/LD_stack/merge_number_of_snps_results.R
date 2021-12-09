out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
library(data.table)
result.list = list()
temp = 1
for(i in 2:5){
  for(l in 1:3){
    for(m in 1:4){
      for(i1 in c(1,4,5,2,3)){
        load(paste0(out.dir,"n_snp_result_eth_",i,"_rho_",l,"_size_",m,"_ga_",i1,".rdata"))        
        result.list[[temp]] = result
        temp = temp + 1
      }
    }
  }
}
library(tidyverse)
result = rbindlist(result.list)
result = result %>% 
  mutate(ga_reorder = case_when(
    ga_vec==1 ~ 1,
    ga_vec==2 ~ 4,
    ga_vec==3 ~ 5,
    ga_vec==4 ~ 3,
    ga_vec==5 ~ 2
  )) %>% 
  select(eth.vec,n.snp.vec,l_vec,m_vec,ga_reorder,ga_vec,method_vec)
result.wide = spread(result,key = "method_vec",value = "n.snp.vec") %>% 
  mutate(cau_vec = case_when(
    l_vec== 1 ~ paste0("0.01"),
    l_vec== 2 ~ paste0("0.001"),
    l_vec== 3 ~ paste0("5E-04")
  ),
  sample_size = case_when(
    m_vec == 1~ "15000",
    m_vec == 2~ "45000",
    m_vec == 3~ "80000",
    m_vec == 4~ "100000"
  ),
  herit_setting = case_when(
    ga_vec==1 ~ "Common SNPs heritatibility as 0.4",
    ga_vec==2 ~ "Whole genome SNPs heritatibility as 0.4",
    ga_vec==3 ~ "Whole genome SNPs heritatibility as 0.4",
    ga_vec==4 ~ "Common SNPs heritatibility as 0.4",
    ga_vec==5 ~ "Common SNPs heritatibility as 0.4",
  ),
  negative_select = case_when(
    ga_vec==1 ~ "Strong negative selection",
    ga_vec==2 ~ "Strong negative selection",
    ga_vec==3 ~ "Strong negative selection",
    ga_vec==4 ~ "No negative selection",
    ga_vec==5 ~ "Mild negative selection",
  ),
  genetic_corre = case_when(
    ga_vec==1 ~ 0.8,
    ga_vec==2 ~ 0.8,
    ga_vec==3 ~ 0.6,
    ga_vec==4 ~ 0.8,
    ga_vec==5 ~ 0.8,
  )) %>% 
  mutate(cau_vec = factor(cau_vec,
                          levels = c("0.01",
                                     "0.001",
                                     "5E-04")),
         sample_size = factor(sample_size,
                              levels = c("15000","45000","80000","100000")))
range(result.wide[,8]/result.wide[,6])
write.csv(result.wide,file = "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/number_of_snps.csv",row.names = F)

