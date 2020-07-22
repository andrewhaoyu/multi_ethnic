#generate temporary sh file to run simulate_LD_data_subset
code <- rep("c",110000)
temp = 1
for(k in 1:1000){
  for(j in 1:22){
    for(i in 1:5){
      if(i==1&j!=2){
        
      }else{
        code[temp] <- paste0("Rscript /data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data_subset.R ",i," ",j," ",k)
        temp = temp+1
      }
      
    }
  }
}
code <- code[1:(temp-1)]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/simulate_LD_data_subset.sh",col.names = F,row.names = F,quote=F)
