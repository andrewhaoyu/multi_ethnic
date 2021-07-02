#goal: preprocess AABC summary level statistics
setwd("/data/zhangh24/multi_ethnic/data/")
library(tidyverse)
library(data.table)
sum.data = fread("./AABC_data/final_metal_4aa_no_ghana1.txt")
colnames(sum.data)[10] = "P"
sum.data = sum.data %>% 
  separate(MarkerName,into = c("CHR","POS","No1","No2"),sep=":",remove=F) %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=F) %>% 
  mutate(ID=MarkerName,
         Effect_allele=toupper(Allele1),
         Alt_allele=toupper(Allele2)) %>% 
select(ID,chr.pos,CHR,POS,Effect_allele,Alt_allele,
         Freq1,FreqSE,Effect,StdErr,P)
  
sum.data.meta = sum.data
l = 1
trait = c("overall","erpos","erneg")
#match SNP in Ghana study
#load ghana
load(paste0("./GBHS_sum/",trait[l],"_sum.rdata"))
sum.data.update = sum.data %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=F) 
#idx <- which(sum.data$POS==114445880)
sum.data.ga = sum.data.update %>% 
  select(chr.pos,Eff_allele,Ref_allele) %>% 
  rename(
         Eff_allele_GA=Eff_allele,
         Ref_allele_GA = Ref_allele) 


sum.data.match = inner_join(sum.data.ga,
                            sum.data.meta,
                            by="chr.pos") %>% 
  filter(((Effect_allele==Eff_allele_GA)&(Alt_allele==Ref_allele_GA))|
           (Effect_allele==Ref_allele_GA)&(Alt_allele==Eff_allele_GA)) 
sum.data.match =sum.data.match %>% 
  mutate(MAF = ifelse(Freq1<=0.5,Freq1,1-Freq1)) %>% 
  select(chr.pos,ID,CHR,POS,Effect_allele,Alt_allele,MAF,Effect,StdErr,P) %>% 
  
sum.data = sum.data.match
sum.data = sum.data %>% 
  mutate(POS = as.numeric(POS),
         CHR=as.numeric(CHR))
  
save(sum.data,file = paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
