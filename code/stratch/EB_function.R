EstimatePrior <- function(beta_tar,sd_tar,
                          beta_eur,sd_eur){
  beta_tar = as.numeric(beta_tar)
  sd_tar = as.numeric(sd_tar)
  beta_eur = as.numeric(beta_eur)
  sd_eur = as.numeric(sd_eur)
  z_tar = beta_tar/sd_tar
  z_eur = beta_eur/sd_eur
  z_mat <-na.omit(cbind(z_tar,z_eur))
  
  
  prior.mat <- cov(z_mat)-diag(2)
  return(prior.mat)
}
EBpost <- function(beta_tar,sd_tar,
                   beta_eur,sd_eur,EBprior){
  
  beta_tar = as.numeric(beta_tar)
  sd_tar = as.numeric(sd_tar)
  beta_eur = as.numeric(beta_eur)
  sd_eur = as.numeric(sd_eur)
  prior.sigma = EBprior
  z_tar = beta_tar/sd_tar
  z_eur = beta_eur/sd_eur
  z_mat <-as.matrix(cbind(z_tar,z_eur))
  sd_mat =  as.matrix(cbind(sd_tar,sd_eur))
  post.sigma = solve(solve(prior.sigma)+diag(2))
  
  z_mat_post = z_mat
  
  p <- ncol(z_mat)
  
  for(k in 1:nrow(z_mat)){
    if(k%%10000==0){print(paste0(k," SNPs completed"))}
    z_temp = z_mat[k,]
    
    #find out nonmissing component
    
    idx <- which(!is.na(z_temp))
    if(length(idx)<p){
      z_temp <- z_temp[idx]
      
      post.sigma_temp = post.sigma[idx,idx,drop=F]
      z_post = post.sigma_temp%*%z_temp
    }else{
      z_post =post.sigma%*%z_temp
    }   
    
    z_mat_post[k,idx] = z_post
  }
  beta_mat_post = z_mat_post
  beta_mat_post[,1] =z_mat_post[,1]*sd_tar
  beta_mat_post[,2] =z_mat_post[,2]*sd_eur
  return(beta_mat_post)
}

