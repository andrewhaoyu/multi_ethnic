method = c("PT","TDLD")
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
for(k in 1:length(method)){
  
  
  system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result; mkdir ",method[k]))
  # for(i in 1:length(eth)){
  #   for(l in 1:length(trait)){
  #     
  #   }
  # }  
}
for(k in 1:length(method)){
  
  
  #system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k]))
  for(i in 1:length(eth)){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k],"; mkdir ",eth[i]))
    # for(l in 1:length(trait)){
    # 
    # }
  }
}


for(k in 1:length(method)){
  
  
  #system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k]))
  for(i in 1:length(eth)){
    #system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k],"; mkdir ",eth[i]))
     for(l in 1:length(trait)){
       system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k],"/",eth[i],"; mkdir ",trait[l]))
     }
  }
}



method = c("PT","TDLD","BestEURPRS_EURcoef",
           "BestEURPRS_tarcoef",
           "BestEURPRS_EBcoef",
           "Weighted_PRS",
           "TDLD_EB")
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
for(k in 1:length(method)){
  
  
  system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/prs; mkdir ",method[k]))
  # for(i in 1:length(eth)){
  #   for(l in 1:length(trait)){
  #     
  #   }
  # }  
}
for(k in 1:length(method)){
  
  
  #system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k]))
  for(i in 1:length(eth)){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/prs/",method[k],"; mkdir ",eth[i]))
    # for(l in 1:length(trait)){
    # 
    # }
  }
}


for(k in 1:length(method)){
  
  
  #system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k]))
  for(i in 1:length(eth)){
    #system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method[k],"; mkdir ",eth[i]))
    for(l in 1:length(trait)){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/prs/",method[k],"/",eth[i],"; mkdir ",trait[l]))
    }
  }
}