library(MetaSubtract)
meta.results = read.table("/Library/Frameworks/R.framework/Versions/4.0/Resources/library/MetaSubtract/extdata/meta_results.txt",header=T)
library(tidyverse)
meta.results.update = meta.results %>% 
  select(chromosome,position,
         reference_allele,other_allele,
         eaf,beta,se,p.value,
         rs_number) %>% 
  mutate(maf = ifelse(eaf<=0.5,eaf,1-eaf))
write.table(meta.results.update,
            "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/MetaSubtract/extdata/meta_results_update.txt",
            row.names = F,
            col.names = T,
            quote=F)
cohort.results = read.table("/Library/Frameworks/R.framework/Versions/4.0/Resources/library/MetaSubtract/extdata/cohort1_results.txt",header=T)
cohort.results.update = cohort.results %>% 
  select(MARKER,CHR,POSITION,EFFECTALLELE,OTHERALLELE,
         BETA,SE,P,EAF)
write.table(cohort.results.update,
          "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/MetaSubtract/extdata/cohort1_results_update.txt",
                                   row.names = F,
                                   col.names = T,
                                   quote=F)
metafile="meta_results_update.txt"
cohortfiles=c("cohort1_results_update.txt")
m1<-meta.subtract(metafile=metafile, cohortfiles=cohortfiles, lambda.meta=1, gc_meta = F,dir=tempdir())
head(m1)
Meta = function(coef_vec,se_vec){
  var_vec = se_vec^2
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  meta_se = sqrt(meta_var)
  return(c(meta_coef,meta_se))
}
coef_vec = c(0.0103790,-0.0086538708)
se_vec = c(0.019930,0.006874788)
Meta(coef_vec,se_vec)




 Metasub = function(coef_meta,se_meta,
                    coef_co,se_co){
   var_meta = se_meta^2
   var_co = se_co^2
   var_co2 = (var_co*var_meta)/(var_co-var_meta)
   coef_co2  = (coef_meta/var_meta-coef_co/var_co)*var_co2
   return(c(coef_co2,sqrt(var_co2)))
 }
 Metasub(-0.006630,0.006499,
         0.0103790,0.019930)
 