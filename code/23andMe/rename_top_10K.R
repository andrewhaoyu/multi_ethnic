
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
library(data.table)
data_dir = "/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/data/23andMe_Top10K/top10K/"
for(i in 1:5){
  for(l in 1:7){
    system(paste0("cd ",data_dir,";",
                  "cd ",eth[i],";",
                  "mv ",trait[l],"_sumdata_top10K.txt ",
                  eth[i],"_",trait[l],"_sumdata_top10K.txt"))    
  }
}
