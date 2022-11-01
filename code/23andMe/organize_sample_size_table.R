#organize sample size table
setwd("/Users/zhangh24/Google Drive (haoyuzhang@hsph.harvard.edu)/23andMe_PRS/Haoyu/cohort_summary")
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")
eth_name = c("European","African American",
             "Latino","East Asian","South Asian")
trait = c("any_cvd","depression",
          "heart_metabolic_disease_burden",
          "height",
          "iqb.sing_back_musical_note",
          "migraine_diagnosis",
          "morning_person")
trait_name = c("Any CVD","Depression",
               "Heart metabolic disease burden",
               "Height",
               "SBMN",
               "Migraine Diagnosis",
               "Morning Person")
total = length(eth_group)*length(trait)
sample.size.table = data.frame(
  eth = rep("c",total),
  trait = rep("c",total),
  train.control = rep(0,total),
  train.case = rep(0,total),
  tun.control = rep(0,total),
  tun.case = rep(0,total),
  vad.control = rep(0,total),
  vad.case = rep(0,total)
)
temp = 1
for(l in 1:length(trait)){
  for(i in 1:length(eth_group)){
    data = read.csv(paste0(eth_group[i],".csv"))
    sample.size.table[temp,1] = eth_name[i]
    sample.size.table[temp,2] = trait_name[l]
    idx <- which(data$pheno==trait[l])
    sample.size.table[temp,3] = data[idx,5]
    sample.size.table[temp,4] = data[idx,4]
    sample.size.table[temp,5] = data[idx,7]
    sample.size.table[temp,6] = data[idx,6]
    sample.size.table[temp,7] = data[idx,3]
    sample.size.table[temp,8] = data[idx,2]
    
    temp = temp + 1
    
  }
}
write.csv(sample.size.table,file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/23andme_sample_size.csv")
sample.size.table$total = rowSums(sample.size.table[,3:8])
library(dplyr)
sample.size.sum = sample.size.table %>% 
  group_by(eth) %>% 
  summarize(average = max(total))
sum(sample.size.sum$average)












sample.size.table %>% filter(trait%in%c("Heart metabolic disease burden",
                                        "Height")==F) %>% 
  group_by(eth) %>% 
  summarize(average = mean(total))

sample.size.table %>% filter(trait%in%c("Heart metabolic disease burden",
                                        "Height")) %>% 
 select(eth,trait,total)



sample.size.table %>% filter(eth!="European") %>% 
  group_by(trait) %>% 
  summarise(per_trait_sum = sum(total)) %>% 
  mutate(mean(per_trait_sum))


sample.size.table %>% filter(eth=="African American"&
                               trait =="Height") 

sample.size.table %>% filter(
                               trait =="Height") 
sample.size.table = sample.size.table %>% 
  mutate(total.cases = train.case + tun.case +  vad.case,
         total.control = train.control + tun.control + vad.control)
sample.size.table %>% filter(
  trait =="Any CVD") 
