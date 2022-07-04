#goal: merge lambda values
eth_name = c("EUR","AFR","HIS","EAS","SAS")
trait <- c("HDL","LDL",
           "logTG",
           "nonHDL",
           "TC")
result = list()
temp = 1
for(l in 1:5){
for(i in 1:5){
  
    load(paste0("/data/zhangh24/multi_ethnic/result/GLGC/lambda_value/lambda_vec_",i,"_",l,".rdata"))
    result[[temp]] = data.frame(Trait = trait[l],
                                eth = eth_name[i],
                                lambda = lambda_vec[1],
                                lambda_1000 = lambda_vec[2])
    temp = temp + 1
  }
}

result.table = rbindlist(result)

write.csv(result.table, file = "/data/zhangh24/multi_ethnic/result/GLGC/lambda_value/lambda_table.csv")
