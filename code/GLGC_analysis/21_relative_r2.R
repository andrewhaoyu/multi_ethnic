#relative r2 performance
me23 = read.csv("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/23andme/prediction_summary.csv")[,-1,drop=F]
me23 = me23 %>% 
  mutate(method_vec = ifelse(method_vec == "PolyPred+", "PolyPred-S+", method_vec))
glgc = read.csv("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/prediction_summary.csv")[,-1,drop=F]
glgc = glgc %>% 
  filter(eth_name!="Latino") %>% 
  select(eth_name,trait_name,method_vec,result)
allofus = read.csv("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/AOU/analysis_result/prediction_summary.csv")[,-1,drop=F]
colnames(me23)[1] = "eth_name"
allofus = allofus %>% 
  mutate(trait_name = ifelse(trait_name == "Height", "Height (AoU)", trait_name)) %>% 
  select(eth_name,trait_name,method_vec,result)
auctosigma2 = function(x){
  return(ifelse(x<=0.5,0, qnorm(x)^2*2))}
result = rbind(me23,glgc,allofus)
me23 = me23 %>% 
  mutate(trait_name = ifelse(trait_name == "Height", "Height (23andMe)", trait_name))

bin_trait =  c(
  "Any cardiovascular disease",
  "Depression",
  "Sing back musical note",
  "Migraine diagnosis",
  "Morning person")
result = rbind(me23,glgc,allofus) %>% 
  mutate(eth_name = ifelse(eth_name=="African American", "African", eth_name)) %>% 
  mutate(result_new = 
           ifelse(trait_name%in%bin_trait,
                  auctosigma2(result),result))
all_trait_name =c("Heart metabolic disease burden",
                  "Height (23andMe)",
                  "Any cardiovascular disease",
                  "Depression",
                  "Sing back musical note",
                  "Migraine diagnosis",
                  "Morning person",
                  "High-density lipoprotein cholesterol",
                  "Low-density lipoprotein cholesterol",
                  "Log triglycerides",
                  "Total cholesterol",
                  "Body mass index",
                  "Height (AoU)")
all_eth_name = c("European", "African", "East Asian", "Latino", "South Asian")
eth_vec = rep(0,100)
trait_vec = rep(0,100)
r2_ratio = rep(0,100)
method_vec = rep(0,100)
temp = 1
r2_tar_vec = rep(0,100)
r2_EUR_vec = rep(0,100)
for(i in 2:5){
  for(l in 1:length(all_trait_name)){
    result_sub_eur = result %>% 
      filter(eth_name == "European" & trait_name == all_trait_name[l])
    #best EUR PRS
    r2_EUR = max(result_sub_eur$result_new)
    
    result_sub_tar = result %>% 
      filter(eth_name == all_eth_name[i] & trait_name == all_trait_name[l])
    #GLGC only have AFR, EAS, SAS; AoU only have AFR
    if(nrow(result_sub_tar)!=0){
      #best target PRS
      r2_tar = max(result_sub_tar$result_new)
      idx <- which.max(result_sub_tar$result_new)
      method_vec[temp] = result_sub_tar$method_vec[idx]
      eth_vec[temp] = all_eth_name[i]
      trait_vec[temp] = all_trait_name[l]
      r2_ratio[temp] = r2_tar/r2_EUR
      r2_tar_vec[temp] = r2_tar
      r2_EUR_vec[temp] = r2_EUR
      temp = temp + 1
      
    }
    
    
    
  }
}
eth_vec = eth_vec[1:(temp-1)]
trait_vec = trait_vec[1:(temp-1)]
method_vec = method_vec[1:(temp-1)]
r2_ratio = r2_ratio[1:(temp-1)]
r2_tar_vec = r2_tar_vec[1:(temp-1)]
r2_EUR_vec = r2_EUR_vec[1:(temp-1)]
ratio_result = data.frame(eth_vec, trait_vec, method_vec, r2_ratio,r2_tar_vec,r2_EUR_vec)
write.csv(ratio_result, file = "/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/relative_r2.csv")
ratio_result %>% group_by(eth_vec) %>% 
  summarize(average = mean(r2_ratio))
