#sample_size
setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/")
data = read.csv("./data/23_sample_size.csv")
library(dplyr)
data = data %>%
  mutate(N_case = ifelse(is.na(as.numeric(N_case)),0,N_case),
         N = as.numeric(N_control)+as.numeric(N_case))
bintrait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

data.sub = data %>% 
  filter(eth=="EUR")
