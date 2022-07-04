#goal: get the sample size of GLGC by ancestry
#UKB is taken out of the samples
library(dplyr)
library(data.table)
data = read.csv("./data/GLGC_raw_sample_size.csv")


data.sub = data %>% 
  filter(Cohort!="UKB")
table(data.sub$Ancestry)

colnames(data.sub)[5] = "N"

data.sub %>% 
  mutate(N = as.numeric(gsub("," , "", N))) %>% 
  group_by(Ancestry) %>% 
  summarise(TotalN = sum(N, na.rm = T))
  


data %>% 
  filter(Ancestry=="Arab living in Kuwait")
