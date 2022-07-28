#goal: merge lambda values
eth_name = c("EUR","AFR","AMR")
trait <- c("height","bmi")
result = list()
temp = 1
for(l in 1:2){
  for(i in 1:3){
    data <- fread(paste0(eth[i],"/",trait[l],".txt"),header=T)
    N.effect = median(data$N)
    load(paste0("/data/zhangh24/multi_ethnic/result/AOU/lambda_value/lambda_vec_",i,"_",l,".rdata"))
    result[[temp]] = data.frame(Trait = trait[l],
                                eth = eth_name[i],
                                N = N.effect,
                                lambda = lambda_vec[1],
                                lambda_1000 = lambda_vec[2])
    temp = temp + 1
  }
}

result.table = rbindlist(result)

write.csv(result.table, file = "/data/zhangh24/multi_ethnic/result/AOU/lambda_value/lambda_table.csv")


