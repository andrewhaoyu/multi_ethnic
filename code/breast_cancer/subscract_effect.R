Metasub = function(coef_meta,se_meta,
                   coef_co,se_co){
  var_meta = se_meta^2
  var_co = se_co^2
  var_co2 = (var_co*var_meta)/(var_co-var_meta)
  coef_co2  = (coef_meta/var_meta-coef_co/var_co)*var_co2
  return(c(coef_co2,sqrt(var_co2)))
}



trait = c("overall","erpos","erneg")
AABC_trait = c("overall","ERpos","ERneg")
for(l in 1:3){
  setwd("/data/zhangh24/multi_ethnic/data/")
  #load ghana
  load(paste0("./GBHS_sum/",trait[l],"_sum.rdata"))
  library(dplyr)
  library(tidyverse)
  sum.data.update = sum.data %>% 
    unite("chr.pos",CHR,POS,sep=":",remove=F) 
  #idx <- which(sum.data$POS==114445880)
  sum.data.ga = sum.data.update %>% 
    select(BETA,SE,MAF,P,chr.pos,Eff_allele,Ref_allele) %>% 
    rename(BETA.GA = BETA,
           SE.GA = SE,
           MAF.GA = MAF,
           P.GA = P) 
  
  
  #load all
  load(paste0("./AABC_data/BC_AFR_",AABC_trait[l],"_meta.rdata"))
  sum.data.meta.update = sum.data.meta %>% 
    unite("chr.pos",CHR,POS,sep=":",remove=F)
  
  sum.data.meta.select = sum.data.meta.update %>% 
    select(ID,chr.pos,CHR,POS,ALT,REF,AF_ALT,beta.ALT,
           SE,P,Rsq_ave,KG.ID) %>% 
    rename(Eff_allele_meta = ALT,
           Ref_allele_meta = REF,
           AF_ALT_meta = AF_ALT,
           BETA_meta = beta.ALT,
           SE_meta = SE,
           P_meta = P)
  
  
  #match snps and alleles
  sum.data.match = inner_join(sum.data.ga,
                              sum.data.meta.select,
                              by="chr.pos") %>% 
    filter(((Eff_allele==Eff_allele_meta)&(Ref_allele==Ref_allele_meta))|
             (Eff_allele==Ref_allele_meta)&(Ref_allele==Eff_allele_meta)) %>% 
    filter(SE.GA>SE_meta)
  sum.data.match.update = sum.data.match %>% 
    mutate(BETA_meta_update = ifelse(Eff_allele==Eff_allele_meta,
                                     BETA_meta,
                                     -BETA_meta))
  BETA_update = sum.data.match.update$BETA_meta_update
  SE_update = sum.data.match.update$SE_meta
  BETA_meta = sum.data.match.update$BETA_meta_update
  SE_meta = sum.data.match.update$SE_meta
  BETA_GA = sum.data.match.update$BETA.GA
  SE_GA = sum.data.match.update$SE.GA
  N = length(BETA_update)
  for(i in 1:N){
    if(i%%1000==0){print(i)}
    temp.result <- Metasub(BETA_meta[i],SE_meta[i],
                           BETA_GA[i],SE_GA[i]
    )
    BETA_update[i] = temp.result[1]
    SE_update[i] = temp.result[2]
  }
  #idx <- which(is.nan(SE_update))
  sum.data.match.update$BETA_update = BETA_update
  sum.data.match.update$SE_update = SE_update
  sum.data.match.update = sum.data.match.update %>% 
    mutate(Z = BETA_update/SE_update,
           P = 2*pnorm(-abs(Z),lower.tail = T),
           MAF = ifelse(AF_ALT_meta<=0.5,AF_ALT_meta,1-AF_ALT_meta)) %>% 
    select(ID,CHR,POS,Eff_allele,Ref_allele,
           Rsq_ave,KG.ID,BETA_update,SE_update,P,MAF
    ) %>% 
    rename(BETA = BETA_update,
           SE = SE_update) 
  sum.data = sum.data.match.update
  save(sum.data,file = paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
  
}
#subscract ghana effect from meta analysis data


# load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
# #setwd("")
# load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.match37_38.rdata")
# mega.chr.pos = fread("/data//zhangh24/multi_ethnic/result/LD_simulation_new/mega_snp_chr_pos.txt",header=F)
