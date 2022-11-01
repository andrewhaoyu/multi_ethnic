setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/")
library(tidyverse)
glgc_sample_size = read.csv("./data/GLGC_sample_size.csv")
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
sum(aou_sample_size_sum$max_size)
aou_sample_size_sum_noEUR = aou_sample_size %>% filter(ethnic != "EUR") %>% 
  group_by(ethnic) %>% summarize(max_size = max(sample_size))
sum(aou_sample_size_sum_noEUR$max_size)
aou_sample_size_update = aou_sample_size %>% rename(eth = ethnic)

ukb_sample_size = read.csv("./data/UKBB_sample_size.csv")
ukb_sample_size = ukb_sample_size %>% 
  filter(ethnic!="EUR") %>% 
  mutate(total = tuning + validation)
ukb_sample_size_sum =  ukb_sample_size %>% group_by(ethnic) %>% 
  summarize(max_size = max(total))
sum(ukb_sample_size_sum$max_size)
ukb_sample_size_sum_noEUR = ukb_sample_size_sum
ukb_sample_size_update = ukb_sample_size %>% select(ethnic, trait, total) %>% 
  rename(sample_size = total, eth = ethnic)

me23_sample_size = read.csv("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/23andme_sample_size.csv")

me23_sample_size = me23_sample_size %>% 
  mutate(total = train.control + train.case+ tun.control+ tun.case+vad.control+vad.case)
me23_sample_size_update = me23_sample_size %>% 
  select(eth,trait, total) %>% 
  mutate(eth_update = case_when(
    eth == "European" ~ "EUR",
    eth == "African American" ~ "AFR",
    eth == "Latino" ~ "AMR",
    eth == "East Asian" ~ "EAS",
    eth == "South Asian" ~ "SAS"
  )) %>% 
  select(eth_update, trait, total) %>% 
  rename(eth = eth_update,
         sample_size = total)
me23_sample_size_sum =  me23_sample_size %>% group_by(eth) %>% 
  summarize(max_size = max(total))
sum(me23_sample_size_sum$max_size)
me23_sample_size_sum_noEUR = me23_sample_size %>% filter(eth != "European") %>% 
  group_by(eth) %>% summarize(max_size = max(total))
sum(me23_sample_size_sum_noEUR$max_size)

sum(glgc_sample_size_sum$max_size)+
  sum(aou_sample_size_sum$max_size)+
  sum(ukb_sample_size_sum$max_size)+
  sum(me23_sample_size_sum$max_size)


sum(glgc_sample_size_sum_noEUR$max_size)+
  sum(aou_sample_size_sum_noEUR$max_size)+
  sum(ukb_sample_size_sum_noEUR$max_size)+ 
  sum(me23_sample_size_sum_noEUR$max_size)



com_data = rbind(glgc_sample_size_update, aou_sample_size_update, ukb_sample_size_update,me23_sample_size_update )
com_data_sum = com_data %>% group_by(eth) %>% 
  summarise(max_N = max(sample_size))





#organize table for glgc and all of us

glgc_sample_size = read.csv("./data/GLGC_sample_size.csv")
eth_vec = c("EUR", "AFR", "AMR", "EAS", "SAS")
eth_name = c("European", "African", "Hispanic", "East Asian", "South Asian")
trait_vec = c("HDL", "LDL", "logTG", "TC")
result_list = list()
temp = 1
for(l in 1:length(trait_vec)){
  for(i in 1:length(eth_vec)){
    ethnic = glgc_sample_size$ethnic
    trait = glgc_sample_size$trait
    sample_size = glgc_sample_size$sample_size
    idx <- which(ethnic==eth_vec[i]& trait==trait_vec[l])
    temp_result = data.frame(trait = trait_vec[l], 
                             eth = eth_name[i], 
                             sample_size = sample_size[idx])
    result_list[[temp]] = temp_result
    temp = temp + 1
  }
}
glgc_table = rbindlist(result_list)
write.csv(glgc_table, file = "./result/GLGC/glgc_sample_size_clean.csv", row.names = F)

aou_sample_size = read.csv("./data/AoU_sample_size.csv")
eth_vec = c("EUR", "AFR", "AMR")
eth_name = c("European", "African", "Hispanic")
trait_name = c("Height", "BMI")
trait_vec = c("height", "bmi")
result_list = list()
temp = 1
for(l in 1:length(trait_vec)){
  for(i in 1:length(eth_vec)){
    ethnic = aou_sample_size$ethnic
    trait = aou_sample_size$trait
    sample_size = aou_sample_size$sample_size
    idx <- which(ethnic==eth_vec[i]& trait==trait_vec[l])
    temp_result = data.frame(trait = trait_name[l], 
                             eth = eth_name[i], 
                             sample_size = sample_size[idx])
    result_list[[temp]] = temp_result
    temp = temp + 1
  }
}
aou_table = rbindlist(result_list)
write.csv(aou_table, file = "./result/AoU/aou_sample_size_clean.csv", row.names = F)

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
