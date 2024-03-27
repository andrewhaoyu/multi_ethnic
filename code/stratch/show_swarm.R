args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

a = i-l
print(a)

# save(a, file = 
#        paste0("/data/zhangh24/a_",i,"_",j,".rdata"))