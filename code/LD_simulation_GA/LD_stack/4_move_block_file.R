#goal: process some of the temporary folders
eth = c("eur", "afr", "amr", "eas", "sas")

for(i in 1:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/block_ld_temp/",eth[i],"/block_ld/; ",
  "rm *.ld; rm -rf snplist_ldblk; mkdir snplist_ldblk; ",
  "cd ../; mv block_ld ../../ldblk_1kg_",eth[i],"/"))  
}


for(i in 1:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",eth[i],"; ",
                "rm *.hdf5"))
}

for(i in 1:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",eth[i],"; ",
                "mv block_ld ../block_ld_temp/",eth[i],"/"))  
}


for(i in 1:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/block_ld_temp/",eth[i],"/ ;",
                "mv block_ld ../../ldblk_1kg_",eth[i],"/"))  
}




for(i in 1:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",eth[i],"/ ;",
                "rm *.hdf5; rm snpinfo_1kg_hm3; cd block_ld; ",
                "rm *ld; cd snplist_ldblk; rm snplist_blk*; rm blk_chr; rm blk_size"))  
}


system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA; ",
              "cp snpinfo_mult_1kg_hm3 /data/zhangh24/test1/ref2/."))

system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA; ",
              "cp snpinfo_mult_1kg_hm3 /data/zhangh24/test1/ref2/."))




for(i in c(1,4)){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"; ",
                "cp ldblk_1kg_chr22.hdf5 /data/zhangh24/test1/ref2/ldblk_1kg_",tolower(eth[i]),"/"))
}


for(i in c(1,4)){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"; ",
         "cp snpinfo_1kg_hm3 /data/zhangh24/test1/ref2/ldblk_1kg_",tolower(eth[i])))
}


for(i in c(1,4)){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_",tolower(eth[i]),"; ",
                "cp snpinfo_1kg_hm3 /data/zhangh24/test1/ref1/ldblk_1kg_",tolower(eth[i])))
}



eth = c("EUR", "AFR", "AMR", "EAS", "SAS")
for(i in 1:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_GA/; ",
                "mkdir ",eth[i],"/prscsx_mega/"))  
}



for(i in c(1,2)){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"; ",
                "cp ldblk_1kg_chr1.hdf5 /data/zhangh24/test1/ref2/ldblk_1kg_",tolower(eth[i]),"/ldblk_1kg_chr1.hdf5"))
}



for(i in 1:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/; ",
                "mkdir xpass"))
}

for(i in 1:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/; ",
                "mkdir polypred"))
  
}



temp = fread("/data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/Polypred/2_run_SBayesR.sh")


idx <- which(temp$V3==2&temp$V4==3&temp$V5==4&temp$V6==1&temp$V7==4)

