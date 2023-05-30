args = commandArgs(trailingOnly = T)
#i1 is the first command line parameter
i1 = as.numeric(args[[1]])

a = i1
print(a)
save(a, file = 
       paste0("/n/holystore01/LABS/xlin/Lab/hzhang/multi_ethnic/result/test/a_",
               i1,".rdata"))
