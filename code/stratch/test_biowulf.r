args = commandArgs(trailingOnly = T)

i1 = as.numeric(args[[1]])

#run linear regression for simulate data
#i1 represent data index (from 1 to 100)
#Y is the outcome
#X is a matrix with 10000 rows and 100 columns

load(paste0("/gpfs/gsfs10/users/BC_risk_prediction/test_biowulf/simulated_data_",i1))
Y = data[[1]]
X = data[[2]]
N = nrow(X)
P = ncol(X)
beta = rep(0,P)
for(l in 1:P){
  model = lm(Y~X[,l])
  beta[l] = coefficients(model)[2]
}
save(beta,file = paste0("/gpfs/gsfs10/users/BC_risk_prediction/test_biowulf/result_",i1))













#simulate data
# N = 10000
# P = 100
# for(i1 in 1:100){
#   print(i1)
#   X = matrix(rnorm(N*P),N,P)
#   Y = 2+X[,1]*0.5
#   
#   data = list(Y,X)
#   save(data,file = paste0("/gpfs/gsfs10/users/BC_risk_prediction/test_biowulf/simulated_data_",i1))
# }
