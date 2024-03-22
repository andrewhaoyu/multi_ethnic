#function to convert OR and P to log-OR (beta) and se

ConvertPtoSE <- function(OR, P){
  #beta is the log-odds ratio
  # Z statistics = beta/se
  # P = 2*pnorm(abs(Z), lower.tail = T)
  #therefore, with beta and p, you can get se
  beta = log(OR)
  Z = qnorm(P/2)
  SE = abs(beta/Z)
  return(c(beta,SE))
}


  
  
#test example
OR = 1.04
se = 0.01
Z = log(OR)/se
P = 2*pnorm(-abs(Z), lower.tail = T)
ConvertPtoSE(OR, P)
  
  
  
