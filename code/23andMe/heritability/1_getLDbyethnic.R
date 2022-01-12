#calculate the LD-score for each ethnic using 1000 genome data
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])

eth = c("EUR","AFR","AMR","EAS","SAS")
# for(i in 1:5){
#   for(j in 1:22){
    sid<-Sys.getenv('SLURM_JOB_ID')
    dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
    temp.dir = paste0('/lscratch/',sid,'/test/',eth[i],"/")
    system(paste0("mkdir ",temp.dir))
    kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
    
    system(paste0("cp ",kg.dir,eth[i],"/chr",j,".bed ",temp.dir,"chr",j,".bed"))
    system(paste0("cp ",kg.dir,eth[i],"/chr",j,".bim ",temp.dir,"chr",j,".bim"))
    system(paste0("cp ",kg.dir,eth[i],"/chr",j,".fam ",temp.dir,"chr",j,".fam"))
    
    #remove HLC region (25002566-34994414)
    library(dplyr)
    if(j==6){
      bim = read.table(paste0(temp.dir,"chr",j,".bim"))
      head(bim)
      bim.hla = bim %>% 
        filter(V4>=25002566&
                 V4<=34994414)
      dim(bim.hla)
      write.table(bim.hla$V2,file = paste0(temp.dir,"remove_region"),row.names = F,col.names = F,quote=F)
      system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,"chr",j," --exclude ",temp.dir,"remove_region --make-bed --out ",temp.dir,"chr_clean",j))
      if(i!=2){
        #AMR 1000 genomes data is relateively limited
        system(paste0("module load ldsc; ", 
                      "cd /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"; ",
                      "python /data/zhangh24/ldsc/ldsc.py ",
                      "--bfile ",temp.dir,"chr_clean",j," ",
                      "--l2 --ld-wind-kb 1000 ",
                      "--out /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/",eth[i],"_ldsc/",j))
        
      }else{
        system(paste0("module load ldsc; ", 
                      "cd /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"; ",
                      "python /data/zhangh24/ldsc/ldsc.py ",
                      "--bfile ",temp.dir,"chr_clean",j," ",
                      "--l2 --ld-wind-kb 2000 ",
                      "--out /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/",eth[i],"_ldsc/",j))
        
      }
      
    }else{
      if(i!=2){
        system(paste0("module load ldsc;", 
                      "cd /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"; ",
                      "python /data/zhangh24/ldsc/ldsc.py ",
                      "--bfile ",temp.dir,"chr",j," ",
                      "--l2 --ld-wind-kb 1000 ",
                      "--out /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/",eth[i],"_ldsc/",j))
      }else{
        system(paste0("module load ldsc;", 
                      "cd /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"; ",
                      "python /data/zhangh24/ldsc/ldsc.py ",
                      "--bfile ",temp.dir,"chr",j," ",
                      "--l2 --ld-wind-kb 2000 ",
                      "--out /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/",eth[i],"_ldsc/",j))
        
      }
      
    }
    system(paste0("rm -rf ",temp.dir))
    #need to create snp list based on the file for ldsc    
#   }
# }
#need to create snp list based on the file for ldsc    
# data.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
# library(data.table)
# for(i in 1:5){
#   data.list = list()
#   for(j in 1:22){
#     bim = fread(paste0(data.dir,eth[i],"/chr",j,".bim"))
#     data.list[[j]] = bim
#   }
#   data = rbindlist(data.list)
#   data.select = data %>%
#     select(V2,V5,V6) %>%
#     rename(SNP=V2,
#            A1 = V5,
#            A2 = V6)
#   fwrite(data.select,file = paste0(data.dir,eth[i],"/snp_list"),row.names = F,col.names = T,quote=F,sep = " ")
# }
# # #create result folder for ldsc
# trait <- c("any_cvd","depression",
#            "heart_metabolic_disease_burden",
#            "height",
#            "iqb.sing_back_musical_note",
#            "migraine_diagnosis",
#            "morning_person")
# result.dir = "/data/zhangh24/multi_ethnic/result/cleaned/herit"
# for(i in 1:5){
#   system(paste0("cd ",result.dir,"; mkdir ",eth[i]))
#   for(l in 1:7){
#     system(paste0("cd ",result.dir,"/",eth[i],"; mkdir ",trait[l]))
#   }
# }
