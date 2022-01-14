load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.rdata")
order.trait =  c("Heart metabolic disease burden",
                                   "Height",
                                   "Any CVD",
                                   "Depression",
                                   "SBMN",
                                   "Migraine Diagnosis",
                                   "Morning Person")
ethname = c("European","African American",
        "Latino","East Asian","South Asian")
methodname = c("C+T",
           "LDpred2",
           "Best EUR SNP (C+T)",
           "Best EUR SNP (LDpred2)",
           "Weighted PRS",
           #                  "PRS-CSx",
           "CT-SLEB (two ethnics)",
           "CT-SLEB (five ethnics)")
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
                 Method==methodname[k]&
                 trait==order.trait[l]) %>% 
        select(eth,trait,Method,result)
      if(nrow(prediction.sub)!=0){
        result.list[[temp]] = prediction.sub
        temp = temp+1
      }
    }
  }
}
order.result = rbindlist(result.list)
write.csv(order.result, file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.csv")
