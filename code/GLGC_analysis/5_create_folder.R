
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

