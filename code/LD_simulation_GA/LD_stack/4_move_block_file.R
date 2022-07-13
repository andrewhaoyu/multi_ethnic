#goal: process some of the temporary folders
eth = c("eur", "afr", "amr", "eas", "sas")

for(i in 2:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/block_ld_temp/",eth[i],"/block_ld/; ",
  "rm *.ld; rm -rf snplist_ldblk; mkdir snplist_ldblk; ",
  "cd ../; mv block_ld ../../ldblk_1kg_",eth[i],"/"))  
}


for(i in 2:5){
  system(paste0("cd /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",eth[i],"; ",
                "rm *.hdf5"))
}