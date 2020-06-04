n = 100
data <- data.frame(id1 = paste0("id",c(1:n)),
                   infor1 = paste0("infor1",c(1:n)),
                   infor2 = paste0("infor2",c(1:n)),
                  x1= rnorm(n),
                   x2 = rnorm(n),
                   x3 = rnorm(n),
                   x4 = rnorm(n),
                   x5 = rnorm(n))
library(dplyr)
var_name = paste0("x",c(1:5))
for(i in 1:5){
  var_name 
  data.new = data %>% filter(get(var_name[i])<=0) %>% 
  #write.table(data.new,file="")
}
