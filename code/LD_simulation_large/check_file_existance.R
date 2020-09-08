eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"




mis_mat <- NULL

total = 4*10*3*4*22*12*12
all.files <- rep("c",total/4)
all_mat <- matrix(0,total/4,7)
#total <- 
for(i in 2:5){
  files <- dir(path = paste0(cur.dir,eth[i],"/prs/"),pattern = "prs_two_dim_(.*).profile")
  print(i)
  temp = 1
  for(i_rep in 1:10){
    print(i_rep)
    for(l in 1:3){
      for(m in 1:4){
        for(j in 1:22){
          for(k1 in 1:12){
            for(k2 in 1:12){
              all.files[temp]  <- paste0("prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".profile")
             ind_vec <- c(i,i_rep,l,m,j,k1,k2)
             all_mat[temp,] <- ind_vec
              temp  = temp+1
            }  
          }
        }
      }
    }
  } 
  idx <- which(all.files%in%files==F)
  mis_mat <- rbind(mis_mat,all_mat[idx,])
  }
  
# mis_mat[temp-1,]
# mis_mat[temp,]
# mis_mat = mis_mat[1:temp,]
save(mis_mat,file ="/data/zhangh24/multi_ethnic/result/LD_simulation_new/mis_mat.rdata")
