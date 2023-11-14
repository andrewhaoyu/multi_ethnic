eth <- c("EUR","AFR","AMR")
eth_name = c("EUR","AFR","AMR")
trait <- c("height","bmi")


system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; mkdir clumping_result"))
system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/clumping_result; mkdir PT"))

for(i in 1:3){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/; mkdir ",eth[i]))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/" ,eth[i],"/; ",
                  "mkdir ",trait[l]))
  }
}

system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; mkdir prs/; cd prs/; mkdir PT"))
for(i in 1:3){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/PT/; mkdir ",eth[i]))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/PT/" ,eth[i],"/; ",
                  "mkdir ",trait[l]))
  }
}

system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; mkdir BestEUR"))
for(i in 2:3){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/BestEUR/; mkdir ",eth[i]))
}

system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/; mkdir BestEUR"))
for(i in 2:3){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/BestEUR/; mkdir ",eth[i]))
}


for(i in 2:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/BestEUR/" ,eth[i],"/; ",
                  "mkdir ",trait[l]))
  }
}

method_vec = c("weighted_prs", "CTSLEB", "XPASS", "PRSCSX", "polypred")

for(k in 1:length(method_vec)){
  method = method_vec[k]
  # system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; mkdir ",method))
  # for(i in 2:3){
  #   system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/",method,"/; mkdir ",eth[i]))
  # }
  # 
  # 
  # for(i in 2:3){
  #   for(l in 1:2){
  #     system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/",method,"/" ,eth[i],"/; ",
  #                   "mkdir ",trait[l]))
  #   }
  # }
  # 
  # system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/; mkdir ",method))
  # for(i in 2:3){
  #   system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/",method,"/; mkdir ",eth[i]))
  # }

  for(i in 2:3){
    for(l in 1:2){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }
  
}





for(k in 1:length(method_vec)){
  method = method_vec[k]
  #system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; mkdir ",method))
  for(i in 1:1){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/",method,"/; mkdir ",eth[i]))
  }


  for(i in 1:1){
    for(l in 1:2){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }

  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/; mkdir ",method))
  for(i in 1:1){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/",method,"/; mkdir ",eth[i]))
  }
  
  for(i in 1:1){
    for(l in 1:2){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }
  
}



system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/herit/"))
for(i in 1:3){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/herit/",eth[i],"/"))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/herit/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}





system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/PT"))
for(i in 1:3){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/herit/",eth[i],"/"))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/herit/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}



system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/CT_SLEB_all"))
for(i in 1:3){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/CT_SLEB_all/",eth[i],"/"))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/CT_SLEB_all/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}





system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSX_all"))
for(i in 1:3){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSX_all/",eth[i],"/"))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSX_all/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}





system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/EURPRS"))
for(i in 1:3){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/EURPRS/",eth[i],"/"))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/EURPRS/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}




system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/weightedprsall"))
for(i in 1:3){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/AOU/prs/weightedprsall/",eth[i],"/"))
}
for(i in 1:3){
  for(l in 1:2){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/prs/weightedprsall/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}





system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/; mkdir boot_result"))




method_vec = c("PT","weighted_prs", "CTSLEB", "PRSCSX", "BESTEUR")

for(k in 1:length(method_vec)){
  method = method_vec[k]
  system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/boot_result/; mkdir ",method))
  
  for(i in 1:3){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/boot_result/",method,"/;",
                  "mkdir ",eth[i]))
    for(l in 1:2){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/AOU/boot_result/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }
  
}



