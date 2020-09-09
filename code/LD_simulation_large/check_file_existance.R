eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"


mis_mat <- NULL

total = 4*10*3*4*22*12*12
all.files <- rep("c",total/4)
all_mat <- matrix(0,total/4,8)
all_mat2 <- matrix(0,total/4,8)

#total <- 
for(i in 2:5){
  files <- dir(path = paste0(cur.dir,eth[i],"/prs/"),pattern = "prs_two_dim_(.*).profile")
  print(i)
  temp = 1
  total = 0
  for(i_rep in 1:10){
    print(i_rep)
    for(j in 1:22){
      load(paste0(cur.dir,"check_mat_rep",i_rep,"_i_",i,"_chr_",j))
      temp2 = nrow(check_mat)
      all_mat[total+(1:temp2),] <- check_mat
      total = total+temp2
    for(l in 1:3){
      for(m in 1:4){
        
          for(k1 in 1:12){
            for(k2 in 1:12){
              all.files[temp]  <- paste0("prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".profile")
             ind_vec <- c(1,i,i_rep,j,l,m,k1,k2)
             all_mat2[temp,] <- ind_vec
              temp  = temp+1
            }  
          }
        }
      }
    }
  } 
  #check consistency
  all.equal(all_mat[,2:8],all_mat2[,2:8])
  idx <- which(all.files%in%files==F)
  colnames(all_mat) <- c("row_bo","i","i_rep","j","l","m","k1","k2")
  mis_mat_temp <- all_mat[idx,]
  #only keep the files with row_bo > 0
  jdx <- which(mis_mat_temp[,1]==1)
  mis_mat <- rbind(mis_mat,mis_mat_temp[jdx,])
  }
  
# mis_mat[temp-1,]
# mis_mat[temp,]
# mis_mat = mis_mat[1:temp,]
save(mis_mat,file ="/data/zhangh24/multi_ethnic/result/LD_simulation_new/mis_mat.rdata")
