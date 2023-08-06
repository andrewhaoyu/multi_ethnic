setwd("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/data/23andMe_Top10K/top10K")
library(openxlsx)

eth <- c("AFR","AMR","EAS","EUR","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

wb <- loadWorkbook("TransPRS_SupplementaryTable_080223_Haoyu.xlsx")
for(l in c(4,5,7)){
for(i in 1:5){
  
    data = read.table(paste0(eth[i],"/",eth[i],"_",trait[l],"_sumdata_top10K.txt"),header=T)
    addWorksheet(wb, paste0(eth[i],"_",trait[l]))
    writeData(wb, sheet = paste0(eth[i],"_",trait[l]), x = data)
  }
}
saveWorkbook(wb, "test_sup.xlsx", overwrite = TRUE)
