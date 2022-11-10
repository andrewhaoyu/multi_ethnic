load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.rdata")
order.trait =  c("Heart metabolic disease burden",
                                   "Height",
                                   "Any CVD",
                                   "Depression",
                                   "SBMN",
                                   "Migraine Diagnosis",
                                   "Morning Person")
trait_name = c("Heart metabolic disease burden",
               "Height",
               "Any cardiovascular disease",
               "Depression",
               "Sing back musical note",
               "Migraine diagnosis",
               "Morning person")
ethname = c("European","African American",
        "Latino","East Asian","South Asian")
methodname = c("CT",
           "LDpred2",
           "Best EUR PRS (CT)",
           "Best EUR PRS (LDpred2)",
           "Weighted PRS (CT)",
           "Weighted PRS (LDpred2)",
           "PolyPred+",
           "XPASS",
           "PRS-CSx",
           "PRS-CSx (five ancestries)",
           "CT-SLEB",
           "CT-SLEB (five ancestries)")
sigma2toauc = function(x){
  ifelse(x==0,0.50,round(pnorm(0.5*sqrt(x)),2))
}
auctosigma2 = function(x){
  return(ifelse(x<=0.5,0, qnorm(x)^2*2))
}
result.list = list()
temp = 1
for(l in 1:length(order.trait)){
  for(i in 1:length(ethname)){
    for(k in 1:length(methodname)){
      prediction.sub = prediction.result %>% 
        filter(eth==ethname[i]&
                 method_vec==methodname[k]&
                 trait==order.trait[l]) %>% 
        mutate(trait_name = trait_name[l]) %>% 
        select(eth,trait_name, method_vec,result)
      if(nrow(prediction.sub)!=0){
        result.list[[temp]] = prediction.sub
        temp = temp+1
      }
    }
  }
}
order.result = rbindlist(result.list)
write.csv(order.result, file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.csv")
