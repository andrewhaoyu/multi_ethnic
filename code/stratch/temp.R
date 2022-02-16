
#we have beta1, se1 for n observations
#we have beta2, se2 for n observations
library(MASS)
library(mvtnorm)
n = 100
beta = rmvnorm(n,mean = c(0,0),sigma = matrix(c(1,0.5,0.5,1),2,2))

beta_hat = matrix(0,n,2)
sigma_data = matrix(c(1/1000,0,0,1/1000),2,2)
for(k in 1:n){
  #sample size as 1000
  
  beta_hat[k,] = rmvnorm(1,mean = beta[k,],sigma = sigma_data)
}

#assume prior for beta as multivariate normal with mean 0 and unknown variance
#get empirical bayes prior
total = matrix(0,2,2)
for(k in 1:n){
  total = total + beta_hat[k,]%*%t(beta_hat[k,])-sigma_data
}
prior_estimate = total/(n-1)

#get posterior mean
post_estimate = matrix(0,n,2)
for(k in 1:n){
  post_estimate[k,] = solve(solve(prior_estimate)+solve(sigma_data))%*%beta_hat[k,]
}

post_estimate1 = post_estimate

#start from z-statistics

z_hat = matrix(0,n,2)
for(k in 1:n){
  z_hat[k,] = beta_hat[k,]%*%solve(sqrt(sigma_data))
}
#assume prior for z as multivariate normal with mean 0 and unknown variance
total = matrix(0,2,2)
for(k in 1:n){
  total = total + z_hat[k,]%*%t(z_hat[k,])-diag(2)
}
prior_estimate = total/(n-1)
#get posterior mean
post_estimate = matrix(0,n,2)
for(k in 1:n){
  post_estimate[k,] = solve(solve(prior_estimate)+solve(diag(2)))%*%z_hat[k,]
}
#transfer back to beta scale
post_estimate2 = post_estimate%*%sqrt(sigma_data)
post_estimate1/post_estimate2