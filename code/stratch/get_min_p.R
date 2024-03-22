#get the code of min p
SNP = paste0("rs",c(1:10))
p_value = matrix(runif(30,0,1),10,3)
data = data.frame(SNP,p_value)
head(data)
p_min = apply(p_value,1,min)
data = data.frame(SNP,p_value,p_min)
