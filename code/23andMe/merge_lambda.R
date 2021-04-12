eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
im <- c("mega_hap","im")
total <- length(eth)*length(trait)*2

lambda_result_vec <- rep(0,total)
lambda_1000_vec <- rep(0,total)
eth_vec <- rep("c",total)
trait_vec <- rep("c",total)
im_vec <- rep("c",total)
temp = 1
for(i3 in  1:length(im)){
for(i1 in 1:length(eth)){
  for(i2 in 1:length(trait)){
    
      load(paste0("/data/zhangh24/multi_ethnic/result/cleaned/lambda_value/lambda_vec_",i1,"_",i2,"_",i3,".rdata"))
      lambda_result_vec[temp] = lambda_vec[1]
      lambda_1000_vec[temp] = lambda_vec[2]
      eth_vec[temp] = eth[i1]
      trait_vec[temp] = trait[i2]
      im_vec[temp] = im[i3]
      temp = temp+1
    }
  }
}

result <- data.frame(eth_vec,trait_vec,
                     im_vec,lambda_result_vec,
                     lambda_1000_vec)

write.csv(result,file = "/data/zhangh24/multi_ethnic/result/cleaned/lambda_value/lambda_table.csv",row.names = F)
