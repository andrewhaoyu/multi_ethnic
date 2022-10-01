n = 1000
z_mat = matrix(rnorm(5*n),n,5)
idx =which(apply(z_mat, 1, function(r) any(abs(r)>=3))==T)
library(dplyr)
test = z_mat %>% filter_all(any_vars(abs(.) >=3))
