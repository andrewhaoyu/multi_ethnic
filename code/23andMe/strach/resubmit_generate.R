eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
phi = c("1e+00","1e-02","1e-04","1e-06")
file.list = list()
temp = 1
for(i in 2:5){
  for(l in 1:7){
    out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/",eth[i],"/",trait[l],"/")
    setwd(out.dir.prs)
    files = dir(out.dir.prs,pattern = "chr")
    for(v in 1:4){
      for(j in 1:22){
        file = paste0("sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt")
        if(file%in%files==F){
          file.list[[temp]] = data.frame(paste0("Rscript /data/zhangh24/multi_ethnic/code/23andMe/PRScsx_run.R ",i," ",l," ",j," ",v))
          temp = temp+1
        } 
        
      }
      
    }
      
  }
}
file.all = rbindlist(file.list)
