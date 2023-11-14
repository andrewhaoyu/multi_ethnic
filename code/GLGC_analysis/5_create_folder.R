
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("HDL","LDL",
           "logTG",
           "TC")

for(i in 5:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/; mkdir ",eth[i]))
}
for(i in 5:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/" ,eth[i],"/; ",
           "mkdir ",trait[l]))
  }
}

system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/; mkdir prs/; cd prs/; mkdir PT"))
for(i in 5:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/PT/; mkdir ",eth[i]))
}
for(i in 5:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/PT/" ,eth[i],"/; ",
                  "mkdir ",trait[l]))
  }
}
system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/; mkdir BestEUR"))
for(i in 2:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/BestEUR/; mkdir ",eth[i]))
}

system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/; mkdir BestEUR"))
for(i in 2:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/BestEUR/; mkdir ",eth[i]))
}


for(i in 2:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/BestEUR/" ,eth[i],"/; ",
                  "mkdir ",trait[l]))
  }
}


method_vec = c("weighted_prs", "CTSLEB", "XPASS", "PRSCSX", "polypred")

for(k in 1:length(method_vec)){
  method = method_vec[k]
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/; mkdir ",method))
  for(i in 1:5){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/",method,"/; mkdir ",eth[i]))
  }


  for(i in 1:5){
    for(l in 1:4){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }

  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/; mkdir ",method))
  for(i in 1:5){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/",method,"/; mkdir ",eth[i]))
  }

  for(i in 1:5){
    for(l in 1:4){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }
  
}


system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/herit/"))
for(i in 1:5){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/herit/",eth[i],"/"))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/herit/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}


system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/PT"))
for(i in 1:5){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/PT/",eth[i],"/"))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/PT/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}



system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/CT_SLEB_all"))
for(i in 1:5){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/CT_SLEB_all/",eth[i],"/"))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/CT_SLEB_all/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}


system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/PRSCSX_all"))
for(i in 1:5){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/PRSCSX_all/",eth[i],"/"))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/PRSCSX_all/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}


system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/EURPRS"))
for(i in 1:5){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/EURPRS/",eth[i],"/"))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/EURPRS/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}



system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/weightedprsall"))
for(i in 1:5){
  system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/prs/weightedprsall/",eth[i],"/"))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/prs/weightedprsall/",eth[i],";",
                  "mkdir ",trait[l],"/"))  
  }
}





eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("HDL","LDL",
           "logTG",
           "TC")



system(paste0("mkdir /data/zhangh24/multi_ethnic/result/GLGC/boot_result/;
cd /data/zhangh24/multi_ethnic/result/GLGC/boot_result/; mkdir PT"))

for(i in 1:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/boot_result/PT/; mkdir ",eth[i]))
}
for(i in 1:5){
  for(l in 1:4){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/boot_result/PT/" ,eth[i],"/; ",
                  "mkdir ",trait[l]))
  }
}


method_vec = c("BestEUR","weighted_prs", "CTSLEB", "PRSCSX")

for(k in 1:length(method_vec)){
  method = method_vec[k]
  system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/boot_result/; mkdir ",method))
  for(i in 1:5){
    system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/boot_result/",method,"/; mkdir ",eth[i]))
  }
  
  
  for(i in 1:5){
    for(l in 1:4){
      system(paste0("cd /data/zhangh24/multi_ethnic/result/GLGC/boot_result/",method,"/" ,eth[i],"/; ",
                    "mkdir ",trait[l]))
    }
  }
}
