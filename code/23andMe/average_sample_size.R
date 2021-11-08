#sample_size
setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/")
data = read.csv("./data/23_sample_size.csv")
library(dplyr)
data = data %>%
  mutate(N_case = ifelse(is.na(as.numeric(N_case)),0,N_case),
         N = (as.numeric(N_control)+as.numeric(N_case))/0.8)
bintrait <- c("any_cvd","depression",
              
              "iqb.sing_back_musical_note",
              "migraine_diagnosis",
              "morning_person")
contriat = c("heart_metabolic_disease_burden",
             "height")

data.sub = data %>% 
  group_by(eth) %>% 
  summarise_at(vars(N),list(name = mean))

  #filter(Disease%in%bintrait) %>% 
data.sub = data %>% filter(Disease%in%contriat) %>% 
  group_by(eth) %>% 
  summarise_at(vars(N),list(name = mean))

data.sub = data %>% filter(Disease%in%bintrait) %>% 
  group_by(eth) %>% 
  summarise_at(vars(N),list(name = mean))
