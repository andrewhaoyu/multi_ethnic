#goal move the simulated genotype data to the shared data directory
eth = c("EUR", "AFR", "AMR", "EAS", "SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <- "/data/BB_Bioinformatics/simulated_multi_ances_genotype_600K/"
for(i in 2:5){
  system(paste0("cd ",cur.dir,eth[i],"; ",
                "mv all_chr.tag.bed ",out.dir,eth[i],"; ",
                "mv all_chr.tag.bim ",out.dir,eth[i],"; ",
                "mv all_chr.tag.fam ",out.dir,eth[i],"; "))
}