args = commandArgs(trailingOnly = T)
#i represent ethnic group
#l represent trait

i = as.numeric(args[[1]])
set.seed(i)
n = 10000
for( i in 1:n){
  rnorm(n)
}