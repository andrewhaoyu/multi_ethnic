trait = c("overall","erpos","erneg")
setwd("/data/zhangh24/multi_ethnic/data/")
for(l in 1:3){
  load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
  library(tidyverse)
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/mega.match.rdata"))
  #load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.match37_38.rdata")
  sum.data.update = sum.data 
  
  mega.match.update = mega.match %>% 
    unite("chr.pos",CHR,POS,sep=":")
  
  mega.sum = inner_join(mega.match.update,
                        sum.data.update,by="chr.pos") 
  mega.sum.update = mega.sum %>% 
    filter(((Effect_allele==REF)&(Alt_allele==ALT))|
             ((Effect_allele==ALT)&(Alt_allele==REF)))
  sum.data = mega.sum.update %>% 
    select(V1,chr.pos,colnames(sum.data))
  save(sum.data,file = paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS_mega.rdata"))
  
}
