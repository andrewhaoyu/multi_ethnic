setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/simresults_jin")
load("R2.ldpred2.mega.GA1.RData")
#load(paste0("LD.clump.result.2DLD.rdata"))
LDpred2.result.list = list()
for(i1 in 1:5){
  load(paste0("R2.ldpred2.mega.GA",i1,".RData"))
  LDpred2.result.list[[i1]] =  R2.ldpred2.mega
}

LDpred2.result.temp = bind_rows(LDpred2.result.list)

LDpred2.result = LDpred2.result.temp %>% 
  mutate(l_vec = rho,
         m_vec = size,
         ga_vec = GA,
         method_vec = "LDpred2",
         eth.vec = race,
         r2.vec = R2) %>% 
  select(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)
save(LDpred2.result,file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/LDpred2.result.rdata")


LDpred.eur.result.list = list()
for(i1 in 1:5){
  load(paste0("R2.ldpred2.EUR.mega.GA",i1,".RData"))
  LDpred.eur.result.list[[i1]] =  R2.ldpred2.EUR.mega
}

LDpredeur.result.temp = bind_rows(LDpred.eur.result.list)

LDpredEUR.result = LDpredeur.result.temp %>% 
  mutate(l_vec = rho,
         m_vec = size,
         ga_vec = GA,
         method_vec = "Best EUR PRS (LDpred2)",
         eth.vec = race,
         r2.vec = R2) %>% 
  select(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)
save(LDpredEUR.result,file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/LDpredEUR.result.rdata")
