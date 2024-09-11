calculate_se <- function(beta, pvalue) {   
  z_abs <- abs(qnorm(pvalue / 2));
  z <- sign(beta) * z_abs ;
  se <- abs(beta / z)   ;
  return(se) }

calculate_se(0,1)
