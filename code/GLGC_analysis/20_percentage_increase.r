library(dplyr)
me23 = read.csv("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.csv")[,-1,drop=F]
glgc = read.csv("/Users/zhangh24/GoogleDrive/multi_ethnic/result/GLGC/analysis_result/prediction_summary.csv")[,-1,drop=F]
glgc = glgc %>% 
  filter(eth_name!="Latino")
allofus = read.csv("/Users/zhangh24/GoogleDrive/multi_ethnic/result/AOU/analysis_result/prediction_summary.csv")[,-1,drop=F]
colnames(me23)[1] = "eth_name"
auctosigma2 = function(x){
     return(ifelse(x<=0.5,0, qnorm(x)^2*2))}
bin_trait =  c(
                           "Any cardiovascular disease",
                           "Depression",
                           "Sing back musical note",
                           "Migraine diagnosis",
                           "Morning person")
result_prscsx = rbind(me23,glgc,allofus) %>% 
  mutate(eth_name = ifelse(eth_name=="African American", "African", eth_name)) %>% 
  filter(method_vec %in% c("PRS-CSx")) %>% 
  mutate(result_new = 
           ifelse(trait_name%in%bin_trait,
                  auctosigma2(result),result))
result_ctsleb = rbind(me23,glgc,allofus) %>% 
  mutate(eth_name = ifelse(eth_name=="African American", "African", eth_name)) %>% 
  filter(method_vec %in% c("CT-SLEB")) %>% 
  mutate(result_new = 
           ifelse(trait_name%in%bin_trait,
                  auctosigma2(result),result))
  
ratio = result_ctsleb$result_new/result_prscsx$result_new
ratio_result = data.frame(result_prscsx[,1:2],ratio = ratio)

ratio_result %>% group_by(eth_name) %>% summarise(average = mean(ratio))


result_prscsx = rbind(me23,glgc,allofus) %>% 
  mutate(eth_name = ifelse(eth_name=="African American", "African", eth_name)) %>% 
  filter(method_vec %in% c("PRS-CSx (five ancestries)", "PRS-CSx (three ancestries)")) %>% 
  mutate(result_new = 
           ifelse(trait_name%in%bin_trait,
                  auctosigma2(result),result))
result_ctsleb = rbind(me23,glgc,allofus) %>% 
  mutate(eth_name = ifelse(eth_name=="African American", "African", eth_name)) %>% 
  filter(method_vec %in% c("CT-SLEB (five ancestries)", "CT-SLEB (three ancestries)")) %>% 
  mutate(result_new = 
           ifelse(trait_name%in%bin_trait,
                  auctosigma2(result),result))

ratio_all = result_ctsleb$result_new/result_prscsx$result_new
ratio_result_all = data.frame(result_prscsx[,1:2],ratio = ratio_all)

ratio_result_all %>% group_by(eth_name) %>% summarise(average = mean(ratio))

ratio_result_summary = rbind(ratio_result, ratio_result_all)
write.csv(ratio_result_summary, file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/GLGC/analysis_result/ct_sleb_prscsx_ratio_summary.csv",
            row.names = F)
