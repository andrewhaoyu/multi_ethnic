load("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/prediction.result.summary.rdata")
order.trait =  c("HDL",
                 "LDL",
                 "logTG",
                 "TC")
library(dplyr)
library(data.table)
trait_name = c("High-density lipoprotein cholesterol",
               "Low-density lipoprotein cholesterol",
               "Log triglycerides",
               "Total cholesterol")
ethname = c("EUR","AFR",
            "AMR","EAS","SAS")
eth_name = c("European", "African", "Latino", "East Asian", "South Asian")
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
order.result = rbindlist(result.list)
order.result = order.result %>% filter(eth_name !="Latino")
write.csv(order.result, file = "/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/prediction_summary.csv")
