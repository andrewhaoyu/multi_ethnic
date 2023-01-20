#merge run time
#TDLD-two ancestries
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/time_mem/"
total = 100
TDLD_SLEB_time_vec = matrix(0,total,4)
for(t_rep in 1:100){
  load(paste0(out.dir,"TDLD_SLEB_trep_",t_rep,".rdata"))  
  TDLD_SLEB_time_vec[t_rep,] =  t(time_vec[,3])
}
time_function = function(x){
  mean(x)/60
}

TDLD_SLEB_time = apply(TDLD_SLEB_time_vec, 2, time_function )
#TDLD-SLEB all ethnic using all ethnic groups data
total = 50
TDLD_SLEBalleth_time_vec = matrix(0,total,4)
for(t_rep in 1:50){
  
  #prs csx split into four sub jobs
  
  
  load(paste0(out.dir,"TDLD_SLEBalleth_trep_",t_rep,".rdata"))
  TDLD_SLEBalleth_time_vec[t_rep,] =  t(time_vec[,3])
  
}
TDLD_SLEBalleth_time = apply(TDLD_SLEBalleth_time_vec, 2, time_function )

#prs-csx two ancestries
prscsx_time_vec = matrix(0,30,4)
for(t_rep in 1:30){
  
  #prs csx split into four sub jobs
  total.time = rep(0,4)
  for(k in 1:4){
    load(paste0(out.dir,"prscsx_new_trep_",t_rep,"_phi_",k,".rdata"))
    total.time = total.time + time_vec[,3]
  }
  
  
  prscsx_time_vec[t_rep,] =  total.time
}
prscsx_time = apply(prscsx_time_vec, 2, time_function )

#prs-csx five ancestries
prscsx_time_vec = matrix(0,30,4)
for(t_rep in 1:30){
  
  #prs csx split into four sub jobs
  total.time = rep(0,4)
  for(k in 1:4){
    load(paste0(out.dir,"prscsx_five_trep_",t_rep,"_phi_",k,".rdata"))
    total.time = total.time + time_vec[,3]
  }
  
  
  prscsx_time_vec[t_rep, ] =  total.time
}
prscsx_time_five = apply(prscsx_time_vec, 2, time_function )


running_time = rbind(TDLD_SLEB_time,TDLD_SLEBalleth_time,prscsx_time,prscsx_time_five)


