out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
load(paste0(out.dir,"LD.clump.result.CT.95CI.rdata"))

ct_result = result.table
ct_result$method_vec = "CT"
ct_result = ct_result %>% 
  select(eth_vec, r2, r2_low, r2_high, l_vec, ga_vec, method_vec)
load(paste0(out.dir,"eur.snp.reult.95CI.rdata"))
library(dplyr)
eur_result = result.table %>% 
  filter(method_vec =="Best EUR SNP (C+T)") %>% 
  mutate(method_vec ="Best EUR SNP (CT)") %>% 
  select(eth_vec,r2.vec,ci_low,ci_high,l_vec,ga_vec,method_vec) %>% 
  rename(r2_low = ci_low,
         r2_high = ci_high,
         r2 =r2.vec)
  load(paste0(out.dir,"weightedprs.CT.95CI.rdata"))
  weighted_prs = result.table %>% 
    mutate(method_vec = "Weighted PRS (CT)" ) %>% 
    rename(r2 = r2.vec,
           r2_low = ci_low,
           r2_high = ci_high) %>% 
    select(eth_vec,r2,r2_low,r2_high,l_vec,ga_vec,method_vec)
load(paste0(out.dir,"prscsx.95CI.rdata"))  
prs_csx = result.table %>% 
  rename(r2 = r2.vec,
         r2_low = ci_low,
         r2_high = ci_high,
         eth_vec = eth.vec) %>% 
  select(eth_vec, r2, r2_low, r2_high,
         l_vec, ga_vec, method_vec)
load( paste0(out.dir,"LD.clump.result.allethtest.EB.95CI.rdata"))
ct_sleb = alleth.EB.result %>% 
  rename(eth_vec = eth.vec,
         r2 = r2.vec,
         r2_low = ci_low,
         r2_high = ci_high) %>% 
  select(eth_vec, r2, r2_low, r2_high,
         l_vec, ga_vec, method_vec)

final_result = rbind(ct_result,
                     eur_result,
                     weighted_prs,
                     prs_csx,
                     ct_sleb)
save(final_result, file = paste0(out.dir,"simulation_result_ci.rdata"))
