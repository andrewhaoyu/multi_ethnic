#merge run time
#TDLD-two ancestries
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/time_mem/"
total = 99
TDLD_SLEB_time_vec = rep(0,total)
for(t_rep in 2:100){
  load(paste0(out.dir,"TDLD_SLEB_trep_",t_rep,".rdata"))  
  TDLD_SLEB_time_vec[t_rep] =  time[3]
}
TDLD_SLEB_time = mean(TDLD_SLEB_time_vec)/60
#TDLD-SLEB all ethnic using all ethnic groups data
TDLD_SLEBalleth_time_vec = rep(0,50)
for(t_rep in 1:50){
  
  #prs csx split into four sub jobs
  
  
  load(paste0(out.dir,"TDLD_SLEBalleth_trep_",t_rep,".rdata"))
  TDLD_SLEBalleth_time_vec[t_rep] = time[3]
  
}
TDLD_SLEBalleth_time = mean(TDLD_SLEBalleth_time_vec)/60

#prs-csx two ancestries
prscsx_time_vec = rep(0,30)
for(t_rep in 1:30){
  
  #prs csx split into four sub jobs
  total.time = 0
  for(k in 1:4){
    load(paste0(out.dir,"prscsx_trep_",t_rep,"_phi_",k,".rdata"))
    total.time = total.time + time[3]
  }
  
  
  prscsx_time_vec[t_rep] =  total.time
}
prscsx_time = mean(prscsx_time_vec)/60

#prs-csx five ancestries
prscsx_time_vec = rep(0,30)
for(t_rep in 1:30){
  
  #prs csx split into four sub jobs
  total.time = 0
  for(k in 1:4){
    load(paste0(out.dir,"prscsx_five_trep_",t_rep,"_phi_",k,".rdata"))
    total.time = total.time + time_list[3]
  }
  
  
  prscsx_time_vec[t_rep] =  total.time
}
prscsx_time_five = mean(prscsx_time_vec)/60


running_time = c(TDLD_SLEB_time,TDLD_SLEBalleth_time,prscsx_time)
