setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/")
library(tidyverse)


#organize table for glgc and all of us

glgc_snp_number= read.csv("./data/GLGC_snp_number.csv")
colnames(glgc_snp_number) = c("ethnic", "trait",
                              "total_snp","subset1","subset2")
eth_vec = c("EUR", "AFR", "AMR", "EAS", "SAS")
eth_name = c("European", "African", "Hispanic", "East Asian", "South Asian")
trait_vec = c("HDL", "LDL", "logTG", "TC")
result_list = list()
temp = 1
for(l in 1:length(trait_vec)){
  for(i in 1:length(eth_vec)){
    ethnic = glgc_snp_number$ethnic
    trait = glgc_snp_number$trait
    total_snp = glgc_snp_number$total_snp
    subset1 = glgc_snp_number$subset1
    subset2 = glgc_snp_number$subset2
    idx <- which(ethnic==eth_vec[i]& trait==trait_vec[l])
    temp_result = data.frame(trait = trait_vec[l], 
                             eth = eth_name[i], 
                             total_snp = total_snp[idx],
                             subset1 = subset1[idx],
                             subset2 = subset2[idx])
    result_list[[temp]] = temp_result
    temp = temp + 1
  }
}
glgc_table = rbindlist(result_list)
write.csv(glgc_table, file = "./result/GLGC/glgc_snp_number_clean.csv", row.names = F)

aou_snp_number = read.csv("./data/AoU_snp_number.csv")
colnames(aou_snp_number) = c("ethnic", "trait",
                              "total_snp","subset1","subset2")
eth_vec = c("EUR", "AFR", "AMR")
eth_name = c("European", "African", "Hispanic")
trait_name = c("Height", "BMI")
trait_vec = c("height", "bmi")
result_list = list()
temp = 1
for(l in 1:length(trait_vec)){
  for(i in 1:length(eth_vec)){
    ethnic = aou_snp_number$ethnic
    trait = aou_snp_number$trait
    total_snp = aou_snp_number$total_snp
    subset1 = aou_snp_number$subset1
    subset2 = aou_snp_number$subset2
    idx <- which(ethnic==eth_vec[i]& trait==trait_vec[l])
    temp_result = data.frame(trait = trait_vec[l], 
                             eth = eth_name[i], 
                             total_snp = total_snp[idx],
                             subset1 = subset1[idx],
                             subset2 = subset2[idx])
    result_list[[temp]] = temp_result
    temp = temp + 1
  }
}
aou_table = rbindlist(result_list)
write.csv(aou_table, file = "./result/AoU/aou_snp_number_clean.csv", row.names = F)

ukb_sample_size = read.csv("./data/UKBB_sample_size.csv")
ukb_sample_size = ukb_sample_size %>% 
  filter(ethnic!="EUR") %>% 
  mutate(total = tuning + validation)
eth_vec = c("EUR", "AFR", "EAS", "SAS")
eth_name = c("European", "African", "East Asian", "South Asian")
trait_name = c("HDL", "LDL", "logTG", "TC", "Height", "BMI")
trait_vec = c("HDL", "LDL", "logTG", "TC", "height", "bmi")
result_list = list()
temp = 1
for(l in 1:4){
  for(i in 2:4){
    ethnic = ukb_sample_size$ethnic
    trait = ukb_sample_size$trait
    tun_sample_size = ukb_sample_size$tuning
    vad_sample_size = ukb_sample_size$validation
    idx <- which(ethnic==eth_vec[i]& trait==trait_vec[l])
    temp_result = data.frame(trait = trait_name[l], 
                             eth = eth_name[i], 
                             tun_sample = tun_sample_size[idx],
                             vad_sample = vad_sample_size[idx])
    result_list[[temp]] = temp_result
    temp = temp + 1
  }
}


for(l in 5:6){
  for(i in 2:2){
    ethnic = ukb_sample_size$ethnic
    trait = ukb_sample_size$trait
    tun_sample_size = ukb_sample_size$tuning
    vad_sample_size = ukb_sample_size$validation
    idx <- which(ethnic==eth_vec[i]& trait==trait_vec[l])
    temp_result = data.frame(trait = trait_name[l], 
                             eth = eth_name[i], 
                             tun_sample = tun_sample_size[idx],
                             vad_sample = vad_sample_size[idx])
    result_list[[temp]] = temp_result
    temp = temp + 1
  }
}
ukb_table =  rbindlist(result_list) 
write.csv(ukb_table, file = "./result/GLGC/ukb_sample_size_clean.csv", row.names = F)
















glgc_sample_size_sum =  glgc_sample_size %>% group_by(ethnic) %>% 
  summarize(max_size = max(sample_size)) 

sum(glgc_sample_size_sum$max_size)
glgc_sample_size_sum_noEUR = glgc_sample_size %>% filter(ethnic != "EUR") %>% 
  group_by(ethnic) %>% summarize(max_size = max(sample_size))
sum(glgc_sample_size_sum_noEUR$max_size)
glgc_sample_size_update = glgc_sample_size %>% rename(eth = ethnic)

aou_sample_size = read.csv("./data/AoU_sample_size.csv")
aou_sample_size_sum =  aou_sample_size %>% group_by(ethnic) %>% 
  summarize(max_size = max(sample_size))
