#Goal: merge prediction single result
#simulation number 100
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')

n.s <- 100
pop.ind <- 1
prediction.result <- matrix(0,n.s,2)
gr <- 2
for(i1 in 1:n.s){
 load(paste0("./multi_ethnic/result/LDP.result_single",i1,"_",pop.ind,"_",gr))
  #load(paste0("./multi_ethnic/result/LDP.result_single",i1,"_",pop.ind))  
  prediction.result[i1,] <- LDP.result
}
colMeans(prediction.result)
