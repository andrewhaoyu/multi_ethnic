eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"




mis_mat <- matrix(0,100000,5)
temp = 1
#total <- 
for(i in 2:5){
  files <- dir(path = paste0(cur.dir,eth[i],"/prs/"),pattern = "prs_two_dim_(.*).profile")
  print(i)
  for(l in 1:3){
    for(m in 1:4){
      for(j in 1:22){
        for(k1 in 1:12){
          for(k2 in 1:12){
            target.file <- paste0("prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".profile")
            if(target.file%in%files==F){
              mis_mat[temp,] = c(i,l,m,k1,k2)
              temp = temp +1
            }
          }  
      }
      }
    }
  }
}
