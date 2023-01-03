method_name = c("EURPRS","PT","weightedprsall","CT_SLEB_all","PRSCSX_all")
for(k in 1:length(method_name)){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs; cp -r ",method_name[k]," ../prs_final"))  
}
system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; zip -r prs_final.zip prs_final"))