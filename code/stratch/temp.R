
#we have beta1, se1 for n observations
#we have beta2, se2 for n observations
library(MASS)
library(mvtnorm)
n = 100
u_vec = rmvnorm(n,mean = c(0,0),sigma = matrix(c(1,0.5,0.5,1),2,2))
maf = matrix(runif(2*n,0.01,0.5),n,2)
beta  = u_vec/sqrt(2*maf*(1-maf))
sigma_vec = (1/(2*maf*(1-maf)))%*%diag(c(0.5,0.1))
#sigma_vec = matrix(0.001,n,2)
beta_hat = matrix(0,n,2)

for(k in 1:n){
  #sample size as 1000
  
  beta_hat[k,] = rmvnorm(1,mean = beta[k,],sigma = diag(sigma_vec[k,]))
}

#assume prior for beta as multivariate normal with mean 0 and unknown variance
#get empirical bayes prior
total = matrix(0,2,2)
for(k in 1:n){
  total = total + beta_hat[k,]%*%t(beta_hat[k,])-diag(sigma_vec[k,])
}
prior_estimate = total/(n-1)

#get posterior mean
post_estimate = matrix(0,n,2)
for(k in 1:n){
  post_estimate[k,] = solve(solve(prior_estimate)+solve(diag(sigma_vec[k,])))%*%solve(diag(sigma_vec[k,]))%*%beta_hat[k,]
}

post_estimate1 = post_estimate

#start from z-statistics
z_hat = matrix(0,n,2)
for(k in 1:n){
  z_hat[k,] = beta_hat[k,]%*%solve(sqrt(diag(sigma_vec[k,])))
}
#assume prior for z as multivariate normal with mean 0 and unknown variance
total = matrix(0,2,2)
for(k in 1:n){
  total = total + z_hat[k,]%*%t(z_hat[k,])-diag(2)
}
prior_estimate = total/(n-1)
prior_estimate_z = prior_estimate
#get posterior mean
post_estimate = matrix(0,n,2)
for(k in 1:n){
  post_estimate[k,] = sqrt(diag(sigma_vec[k,]))%*%prior_estimate%*%solve(prior_estimate+diag(2))%*%z_hat[k,]
  
}
#transfer back to beta scale
post_estimate2 = post_estimate
post_estimate2/post_estimate1
#post_estimate2 = post_estimate%*%sqrt(diag(sigma_vec[k,]))






#start from u-statistics
#suppose we have allele frequencies in two populations
u_hat = beta_hat*sqrt(2*maf*(1-maf))
sigma_u_vec = sigma_vec*(2*maf*(1-maf))

#assume prior for beta as multivariate normal with mean 0 and unknown variance
#get empirical bayes prior
total = matrix(0,2,2)
for(k in 1:n){
  total = total + u_hat[k,]%*%t(u_hat[k,])-diag(sigma_u_vec[k,])
}
prior_estimate = total/(n-1)
prior_estimate_u = prior_estimate
#get posterior mean
post_estimate = matrix(0,n,2)
for(k in 1:n){
  post_estimate[k,] = solve(solve(prior_estimate)+solve(diag(sigma_u_vec[k,])))%*%solve(diag(sigma_u_vec[k,]))%*%u_hat[k,]
}
post_estimate3 = post_estimate/sqrt(2*maf*(1-maf))

post_estimate1/post_estimate3
post_estimate2/post_estimate3


sigma_mat = matrix(c(0.5,0,0,0.1),2,2)
sigma_half = sqrt(sigma_mat)

sigma_half%*%prior_estimate_z%*%sigma_half
prior_estimate_u


æsolve(sigma_half)%*%u_hat[1,]



u_hat[k,]%*%t(u_hat[k,])-diag(sigma_u_vec[k,])

total = matrix(0,2,2)
for(k in 1:n){
  total = total +sigma_half%*%(z_hat[k,]%*%t(z_hat[k,])-diag(2))%*%sigma_half
}
print(total)
total/(n-1)

              