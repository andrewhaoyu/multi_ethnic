#prepare genotype data
# eth = c("EUR", "AFR", "AMR", "EAS", "SAS")
# 
# system("cp -r /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/ /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/UKBB/.")
# system("cd /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/UKBB/; mv all_data genotype;")
# 
# #prepare phenotype data
# system("cp -r /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/UKBB/.")
# 
# #compress
# system("cd /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/; zip -r UKBB.zip UKBB")
# 
# 
# #prepare GLGC summary statistics
# for(i in 1:5){
#   system(paste0("cp -r /dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/", eth[i], "/ /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/GLGC_cleaned/."))
# }
# 
# system(paste0("cp -r /dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/ /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/GLGC_cleaned/."))
# 
# system(paste0("cd /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/; zip -r GLGC_cleaned.zip GLGC_cleaned"))


#copy the results to prs_final
#compress the file
method_name = c("EURPRS","PT","weightedprsall","CT_SLEB_all","PRSCSX_all")
for(k in 1:length(method_name)){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs; cp -r ",method_name[k]," ../prs_final"))  
}
system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/; zip -r prs_final.zip prs_final"))
