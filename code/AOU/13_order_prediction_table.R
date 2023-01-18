load("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/AOU/analysis_result/prediction.result.summary.rdata")
order.trait =  c("bmi",
                 "height")
trait_name = c("Body mass index",
               "Height")
ethname = c("EUR","AFR",
            "AMR")
eth_name = c("European", "African", "Latino")
methodname = c("CT",
               "LDpred2",
               "Best EUR PRS (CT)",
               "Best EUR PRS (LDpred2)",
               "Weighted PRS (CT)",
               "Weighted PRS (LDpred2)",
               "PolyPred+",
               "XPASS",
               "PRS-CSx",
               "PRS-CSx (three ancestries)",
               "CT-SLEB",
               "CT-SLEB (three ancestries)")
sigma2toauc = function(x){
  ifelse(x==0,0.50,round(pnorm(0.5*sqrt(x)),2))
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
        mutate(trait_name = trait_name[l],
               eth_name = eth_name[i]) %>% 
        select(eth_name,trait_name, method_vec,result, r2_low, r2_high)
      if(nrow(prediction.sub)!=0){
        result.list[[temp]] = prediction.sub
        temp = temp+1
      }
    }
  }
}
order.result = rbindlist(result.list) %>% filter(eth_name != "Latino") 

write.csv(order.result, file = "/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/AOU/analysis_result/prediction_summary.csv")
