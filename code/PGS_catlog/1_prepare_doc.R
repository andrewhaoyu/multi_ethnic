#prepare documentation for PGS catalog
#score file
eth_vec <- c("EUR","AFR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
            "PRSCSx","CTSLEB")
trait_name = c("High-density lipoprotein cholesterol",
                            "Low-density lipoprotein cholesterol",
                            "Log triglycerides",
                            "Total cholesterol")

Details = c("clumping window size=500kb, clumping r2=0.1, p-value threshold = c(5E-08,5E-07,...,0.5,1)",
            "causal SNPs protions vary with length 17 that are evened spaced on a logorithmic sale from 1E-04 to 1",
            "Linear weight combination of best single ancestry LDpred2",
            "gamma-gamma prior hyperparameters a and b to default values of 1 and 0.5, respectively. Parameter phi is varied over the default set of values 1E-06, 1E-04), 1E-02), and 1",
            "w_b= (50kb, 100kb), , clumping r2=(0.01, 0.05, 0.1, 0.2, 0.5, 0.8) , clumping window size = w_b/r2, p-value threshold = c(5E-08,5E-07,...,0.5,1),
            super leaner method was set to be glmnet and ridge to obtain linear weight")

ID = rep("c",100000)
Reported_trait = rep("c",100000)
Score_Method = rep("c",100000)
Score_Details = rep("c",100000)
Number_of_var = rep(0,100000)
GLGC_sample_size = fread("/data/zhangh24/multi_ethnic/data/GLGC_sample_size.csv")
pgs_path = "/data/zhangh24/multi_ethnic/result/GLGC/pgs_catalog"
line_count <- function(pgs_path, id_name) {
  cmd <- paste0("zcat ", file.path(pgs_path, paste0(id_name, ".txt.gz")), " | wc -l")
  return(as.numeric(system(cmd, intern = TRUE))-1)
}

temp = 1
for(l in 1:4){
  for(i in 1:4){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }else if(i==1 & k<=2){
        #EUR only have CT and LDpred2
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }
      
    }
  }
}
#AOU
eth_vec <- c("EUR","AFR")
trait_vec <-c("height","bmi")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
trait_name = c("Height",
               "Body mass index")

for(l in 1:2){
  for(i in 1:2){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }else if(i==1 & k<=2){
        #EUR only have CT and LDpred2
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }
      
    }
  }
}

final = temp - 1
ID = ID[1:final]
Reported_trait = Reported_trait[1:final]
Score_Details = Score_Details[1:final]
Score_Method = Score_Method[1:final]
Number_of_var = Number_of_var[1:final]

PGS_score_data = data.frame(ID,Reported_trait,Score_Method,Score_Details,
                            Number_of_var)
write.csv(PGS_score_data, file = "/data/zhangh24/multi_ethnic/result/GLGC/pgs_catlog_infor/PGS_score_data.csv")





#sample description
eth_vec <- c("EUR","AFR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
trait_name = c("High-density lipoprotein cholesterol",
               "Low-density lipoprotein cholesterol",
               "Log triglycerides",
               "Total cholesterol")
Ancestry_vec = c("European",
                       "African American",
                       "East Asian",
                       "South Asian")
multi_ans = "European, African American, Latino, East Asian, South Asian"
Age_vec = c("Mean = 56.86 years (SD = 8.05)",
            "Mean = 51.82 years (SD = 8.06)",
            "Mean = 52.47 years (SD =7.81)",
            "Mean = 53.37 (SD = 8.43)")

Sex_vec = c(0.47,0.42,0.32,0.53)

ID = rep(NA,100000)
Study_Stage = rep(NA,100000)
PubMed = rep(NA,100000)
N_sample = rep(NA,10000)
Broad_ancestry_infor = rep(NA,10000)
Ancestry_infor = rep(NA,10000)
age_infor = rep(NA,10000)
sex_infor = rep(NA,10000)
GLGC_sample_size = fread("/data/zhangh24/multi_ethnic/data/GLGC_sample_size.csv")
UKB_sample_size = fread("/data/zhangh24/multi_ethnic/data/UKBB_sample_size.csv")
AOU_sample_size = fread("/data/zhangh24/multi_ethnic/data/AoU_sample_size.csv")
temp = 1
for(l in 1:4){
  for(i in 1:4){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        #Variant Association
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Variant association"
        PubMed[temp] = 34887591
        
        if(k%in%c(1,2)){
          Broad_ancestry_infor[temp] = Ancestry_vec[i]
          idx = which(GLGC_sample_size$ethnic == eth & 
                        GLGC_sample_size$trait == trait)
          N_sample[temp] = GLGC_sample_size$sample_size[idx]
        }else{
          Broad_ancestry_infor[temp] = "Multiple Ancestries (report in Additional Ancestry Description column)"
          Ancestry_infor[temp] = multi_ans
          idx = which(GLGC_sample_size$trait == trait)
          N_sample[temp] = sum(GLGC_sample_size$sample_size[idx])
        }
        temp = temp + 1  
        #Score development
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Score development"
        Broad_ancestry_infor[temp] = Ancestry_vec[i]
        idx = which(UKB_sample_size$ethnic == eth & 
                      UKB_sample_size$trait == trait)
        N_sample[temp] = UKB_sample_size$tuning[idx]
        age_infor[temp] = Age_vec[i]
        sex_infor[temp] = Sex_vec[i]
        temp = temp + 1  
        #Testing
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Testing"
        Broad_ancestry_infor[temp] = Ancestry_vec[i]
        idx = which(UKB_sample_size$ethnic == eth & 
                      UKB_sample_size$trait == trait)
        N_sample[temp] = UKB_sample_size$tuning[idx]
        age_infor[temp] = Age_vec[i]
        sex_infor[temp] = Sex_vec[i]
        temp = temp + 1  
        
        }else if(i==1 & k<=2){
          #Variant Association
          id_name = paste0(trait,"_",eth,"_",method)
          ID[temp] = id_name
          Study_Stage[temp] = "Variant association"
          PubMed[temp] = 34887591
          idx = which(GLGC_sample_size$ethnic == eth & 
                        GLGC_sample_size$trait == trait)
          N_sample[temp] = GLGC_sample_size$sample_size[idx]
          if(k%in%c(1,2)){
            Broad_ancestry_infor[temp] = Ancestry_vec[i]
          }else{
            Broad_ancestry_infor[temp] = Ancestry_vec[i]
          }
          temp = temp + 1  
          #Score development
          id_name = paste0(trait,"_",eth,"_",method)
          ID[temp] = id_name
          Study_Stage[temp] = "Score development"
          Broad_ancestry_infor[temp] = Ancestry_vec[i]
          idx = which(UKB_sample_size$ethnic == eth & 
                        UKB_sample_size$trait == trait)
          N_sample[temp] = UKB_sample_size$tuning[idx]
          age_infor[temp] = Age_vec[i]
          sex_infor[temp] = Sex_vec[i]
          temp = temp + 1  
          #Testing
          id_name = paste0(trait,"_",eth,"_",method)
          ID[temp] = id_name
          Study_Stage[temp] = "Testing"
          Broad_ancestry_infor[temp] = Ancestry_vec[i]
          idx = which(UKB_sample_size$ethnic == eth & 
                        UKB_sample_size$trait == trait)
          N_sample[temp] = UKB_sample_size$tuning[idx]
          age_infor[temp] = Age_vec[i]
          sex_infor[temp] = Sex_vec[i]
          temp = temp + 1  
      }
      
    }
  }
}
#AOU
eth_vec <- c("EUR","AFR")
trait_vec <-c("height","bmi")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
Ancestry_vec = c("European",
                 "African American")

Age_vec = c("Mean = 56.86 years (SD = 8.05)",
            "Mean = 51.82 years (SD = 8.06)",
            "Mean = 52.47 years (SD =7.81)",
            "Mean = 53.37 (SD = 8.43)")


Sex_vec_aou = c(0.39,0.42,0.40)
Age_vec_aou = c("Mean = 57.97 years (SD = 17.39)",
            "Mean = 55.43 years (SD = 14.81)",
            "Mean = 46.67 years (SD = 16.97)")
multi_ans = c("European, African American, Latino")
for(l in 1:2){
  for(i in 1:2){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        #Variant Association
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Variant association"
        idx = which(AOU_sample_size$ethnic == eth & 
                      AOU_sample_size$trait == trait)
        N_sample[temp] = AOU_sample_size$sample_size[idx]
        if(k%in%c(1,2)){
          Broad_ancestry_infor[temp] = Ancestry_vec[i]
          idx = which(AOU_sample_size$ethnic == eth & 
                        AOU_sample_size$trait == trait)
          N_sample[temp] = AOU_sample_size$sample_size[idx]
        }else{
          Broad_ancestry_infor[temp] = "Multiple Ancestries (report in Additional Ancestry Description column)"
          Ancestry_infor[temp] = multi_ans
          idx = which(AOU_sample_size$trait == trait)
          N_sample[temp] = sum(AOU_sample_size$sample_size[idx])
        }
        age_infor[temp] = Age_vec_aou[i]
        sex_infor[temp] = Sex_vec_aou[i]
        temp = temp + 1  
        #Score development
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Score development"
        Broad_ancestry_infor[temp] = Ancestry_vec[i]
        idx = which(UKB_sample_size$ethnic == eth & 
                      UKB_sample_size$trait == trait)
        N_sample[temp] = UKB_sample_size$tuning[idx]
        age_infor[temp] = Age_vec[i]
        sex_infor[temp] = Sex_vec[i]
        temp = temp + 1  
        #Testing
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Testing"
        Broad_ancestry_infor[temp] = Ancestry_vec[i]
        idx = which(UKB_sample_size$ethnic == eth & 
                      UKB_sample_size$trait == trait)
        N_sample[temp] = UKB_sample_size$tuning[idx]
        age_infor[temp] = Age_vec[i]
        sex_infor[temp] = Sex_vec[i]
        temp = temp + 1  
        
      }else if(i==1 & k<=2){
        #Variant Association
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Variant association"
        idx = which(AOU_sample_size$ethnic == eth & 
                      AOU_sample_size$trait == trait)
        N_sample[temp] = AOU_sample_size$sample_size[idx]
        
        if(k%in%c(1,2)){
          Broad_ancestry_infor[temp] = Ancestry_vec[i]
        }else{
          Broad_ancestry_infor[temp] = Ancestry_vec[i]
        }
        age_infor[temp] = Age_vec_aou[i]
        sex_infor[temp] = Sex_vec_aou[i]
        temp = temp + 1  
        #Score development
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Score development"
        Broad_ancestry_infor[temp] = Ancestry_vec[i]
        idx = which(UKB_sample_size$ethnic == eth & 
                      UKB_sample_size$trait == trait)
        N_sample[temp] = UKB_sample_size$tuning[idx]
        age_infor[temp] = Age_vec[i]
        sex_infor[temp] = Sex_vec[i]
        temp = temp + 1  
        #Testing
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Study_Stage[temp] = "Testing"
        Broad_ancestry_infor[temp] = Ancestry_vec[i]
        idx = which(UKB_sample_size$ethnic == eth & 
                      UKB_sample_size$trait == trait)
        N_sample[temp] = UKB_sample_size$tuning[idx]
        age_infor[temp] = Age_vec[i]
        sex_infor[temp] = Sex_vec[i]
        temp = temp + 1  
      }
      
      
    }
  }
}

final = temp - 1
ID = ID[1:final]
Study_Stage = Study_Stage[1:final]
PubMed = PubMed[1:final]
N_sample = N_sample[1:final]
age_infor = age_infor[1:final]
sex_infor = sex_infor[1:final]
Ancestry_infor = Ancestry_infor[1:final]
Broad_ancestry_infor = Broad_ancestry_infor[1:final]
PGS_sample_data = data.frame(ID,Study_Stage,PubMed,N_sample,
                            age_infor,sex_infor,Ancestry_infor,Broad_ancestry_infor)
write.csv(PGS_sample_data, file = "/data/zhangh24/multi_ethnic/result/GLGC/pgs_catlog_infor/PGS_sample_data.csv",na = "")







#performance description
prepare documentation for PGS catalog
#score file
eth_vec <- c("EUR","AFR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
trait_name = c("High-density lipoprotein cholesterol",
               "Low-density lipoprotein cholesterol",
               "Log triglycerides",
               "Total cholesterol")

Details = c("clumping window size=500kb, clumping r2=0.1, p-value threshold = c(5E-08,5E-07,...,0.5,1)",
            "causal SNPs protions vary with length 17 that are evened spaced on a logorithmic sale from 1E-04 to 1",
            "Linear weight combination of best single ancestry LDpred2",
            "gamma-gamma prior hyperparameters a and b to default values of 1 and 0.5, respectively. Parameter phi is varied over the default set of values 1E-06, 1E-04), 1E-02), and 1",
            "w_b= (50kb, 100kb), , clumping r2=(0.01, 0.05, 0.1, 0.2, 0.5, 0.8) , clumping window size = w_b/r2, p-value threshold = c(5E-08,5E-07,...,0.5,1),
            super leaner method was set to be glmnet and ridge to obtain linear weight")

ID = rep("c",100000)
Reported_trait = rep("c",100000)
Score_Method = rep("c",100000)
Score_Details = rep("c",100000)
Number_of_var = rep(0,100000)
GLGC_sample_size = fread("/data/zhangh24/multi_ethnic/data/GLGC_sample_size.csv")
pgs_path = "/data/zhangh24/multi_ethnic/result/GLGC/pgs_catalog"
line_count <- function(pgs_path, id_name) {
  cmd <- paste0("zcat ", file.path(pgs_path, paste0(id_name, ".txt.gz")), " | wc -l")
  return(as.numeric(system(cmd, intern = TRUE))-1)
}

temp = 1
for(l in 1:4){
  for(i in 1:4){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }else if(i==1 & k<=2){
        #EUR only have CT and LDpred2
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }
      
    }
  }
}
#AOU
eth_vec <- c("EUR","AFR")
trait_vec <-c("height","bmi")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
trait_name = c("Height",
               "Body mass index")

for(l in 1:2){
  for(i in 1:2){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }else if(i==1 & k<=2){
        #EUR only have CT and LDpred2
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        Reported_trait[temp] = trait_name[l]
        Score_Method[temp] = method
        Score_Details[temp] = Details[k]
        Number_of_var[temp] = line_count(pgs_path, id_name)
        temp = temp + 1  
      }
      
    }
  }
}

final = temp - 1
ID = ID[1:final]
Reported_trait = Reported_trait[1:final]
Score_Details = Score_Details[1:final]
Score_Method = Score_Method[1:final]
Number_of_var = Number_of_var[1:final]

PGS_score_data = data.frame(ID,Reported_trait,Score_Method,Score_Details,
                            Number_of_var)
write.csv(PGS_score_data, file = "/data/zhangh24/multi_ethnic/result/GLGC/pgs_catlog_infor/PGS_score_data.csv")




#performance matrix
eth_vec <- c("EUR","AFR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
method_update_vec = c("CT","LDpred2","Weighted PRS (LDpred2)",
                      "PRS-CSx (five ancestries)",
                      "CT-SLEB (five ancestries)")
trait_name = c("High-density lipoprotein cholesterol",
               "Low-density lipoprotein cholesterol",
               "Log triglycerides",
               "Total cholesterol")


ID = rep("c",100000)
R2 = rep("c",100000)
load("/data/zhangh24/multi_ethnic/result/GLGC/prediction.result.summary.rdata")

temp = 1
for(l in 1:4){
  for(i in 1:4){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        idx <- which(prediction.result$eth==eth&
                       prediction.result$trait==trait&
                       prediction.result$method==method_update_vec[k])
        R2[temp] = prediction.result$result[idx]
       temp = temp + 1  
      }else if(i==1 & k<=2){
        #EUR only have CT and LDpred2
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        idx <- which(prediction.result$eth==eth&
                       prediction.result$trait==trait&
                       prediction.result$method==method_update_vec[k])
        R2[temp] = prediction.result$result[idx]
        temp = temp + 1  
      }
      
    }
  }
}
#AOU
eth_vec <- c("EUR","AFR")
trait_vec <-c("height","bmi")
method_vec <- c("CT","LDpred2","weighted_LDpred2",
                "PRSCSx","CTSLEB")
trait_name = c("Height",
               "Body mass index")
method_update_vec = c("CT","LDpred2","Weighted PRS (LDpred2)",
                      "PRS-CSx (three ancestries)",
                      "CT-SLEB (three ancestries)")
load("/data/zhangh24/multi_ethnic/result/AOU/prediction.result.summary.rdata")
for(l in 1:2){
  for(i in 1:2){
    for(k in 1:length(method_vec)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      method = method_vec[k]
      #EUR only have CT and LDpred2
      if(i!=1){
        ID[temp] = id_name
        idx <- which(prediction.result$eth==eth&
                       prediction.result$trait==trait&
                       prediction.result$method==method_update_vec[k])
        R2[temp] = prediction.result$result[idx]
        temp = temp + 1  
      }else if(i==1 & k<=2){
        #EUR only have CT and LDpred2
        id_name = paste0(trait,"_",eth,"_",method)
        ID[temp] = id_name
        idx <- which(prediction.result$eth==eth&
                       prediction.result$trait==trait&
                       prediction.result$method==method_update_vec[k])
        R2[temp] = prediction.result$result[idx]
        temp = temp + 1  
      }
      
    }
  }
}

final = temp - 1
ID = ID[1:final]
R2 = R2[1:final]
PGS_score_data = data.frame(ID,R2)
write.csv(PGS_score_data, file = "/data/zhangh24/multi_ethnic/result/GLGC/pgs_catlog_infor/PGS_performance_data.csv")
