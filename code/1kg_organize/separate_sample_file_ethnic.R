#Goal: separate sample files of 1000 genome by ethnic group
#so that we can use shapeit to split the haplotype file by ethnic group

sample <- read.table("/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3.sample",header=T)
country <- names(table(sample$GROUP))
library(dplyr)
for(i in 1:length(country)){
  sample.sub = sample %>% 
    filter(GROUP==country[i])
  write.table(sample.sub, file = 
                paste0("/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/",country[i],".sample"),
              row.names = F,
              col.names = T,
              quote = F)
  
}

#Goal: unzip the 1000 genome data
a <- rep("c",44)
temp <- 1
for(i in 1:22){
  a[temp] <- paste0("gunzip /spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".hap.gz")
  temp <- temp+1
}

for(i in 1:22){
  a[temp] <- paste0("gunzip /spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend.gz")
  temp <- temp+1
}

write.table(a,file = "/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/gunzip_all_hap_lengend.sh",
            row.names=F,
            col.names=F,
            quote=F)



