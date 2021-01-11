pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres),stringsAsFactors = F)
temp = 1

for(k in 1:length(pthres)){
  
  #  select(SNP,A1,BETA)
  
  
  if(nrow(prs.file)>0){
    q_range[temp,1] = paste0("p_value_",k)
    q_range[temp,3] = pthres[k]
    temp = temp+1
  }
  
  
}
q_range = q_range[1:(temp-1),]

for(i in 1:2){
  for(l in 1:3){
    
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")

write.table(q_range,file = paste0(temp.dir,"q_range_file"),row.names = F,col.names = F,quote=F)
  }
}