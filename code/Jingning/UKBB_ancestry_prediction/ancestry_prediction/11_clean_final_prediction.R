
library(readr)
library(stringr)

final <- read_tsv("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/predicted_ancestry_all.txt")

m <- final$predict_ancestry == "EUR"; m[is.na(m)] <- F
final_nonEUR <- final[!m,]

for(ethnic in c("AFR","AMR","EAS","SAS")){
  tmp <- format(final_nonEUR$IID[final_nonEUR$predict_ancestry == ethnic], trim = T, scientific = F)
  print(length(tmp))
  write_tsv(data.frame(FID=tmp, IID=tmp), paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/",ethnic,".id"), col_name=F)
}


