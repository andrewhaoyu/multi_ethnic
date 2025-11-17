eth = c("EUR", "AFR", "AMR", "EAS", "SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
for(i in 1:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"; ",
                "cp clump_ref_all_chr.bed /data/zhangh24/multi_ethnic/ref_genotype/",eth[i],"/. ;",
                "cp clump_ref_all_chr.bim /data/zhangh24/multi_ethnic/ref_genotype/",eth[i],"/. ;",
                "cp clump_ref_all_chr.fam /data/zhangh24/multi_ethnic/ref_genotype/",eth[i],"/. ;",
                "cp all_chr_test.mega.bed /data/zhangh24/multi_ethnic/ref_genotype/",eth[i],"/. ;",
                "cp all_chr_test.mega.bim /data/zhangh24/multi_ethnic/ref_genotype/",eth[i],"/. ;",
                "cp all_chr_test.mega.fam /data/zhangh24/multi_ethnic/ref_genotype/",eth[i],"/. ;"))
}
system(paste0("cd /data/zhangh24/multi_ethnic/; ",
              "zip -r ref_genotype.zip ref_genotype"))