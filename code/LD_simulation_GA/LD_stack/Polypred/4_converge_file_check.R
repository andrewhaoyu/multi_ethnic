#check whether the sbayesR converged
mis_vec_list = list()

mis_temp =  1

for(i in 1:1){
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
  
  for(l in 1:3){
    for(m in 4:4){
      for(i_rep in 1:10){
        for(i1 in 1:5){
          files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes"),full.names = T)
          out_file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes")
          if(out_file%in%files==F){
            mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
            mis_temp = mis_temp + 1
          }else{
            temp_str = system(paste0("wc -l ",out_file), intern= T)
            num = as.integer(gsub(out_file, "", temp_str))
            if(num == 0){
              mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
              mis_temp = mis_temp + 1
            }
          }
          
        }
      }
    }
  }
}

for(i in 2:5){
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
  
  for(l in 1:3){
    for(m in 1:4){
      for(i_rep in 1:10){
        for(i1 in 1:5){
          files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes"),full.names = T)
          out_file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes")
          if(out_file%in%files==F){
            mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
            mis_temp = mis_temp + 1
          }else{
            temp_str = system(paste0("wc -l ",out_file), intern= T)
            num = as.integer(gsub(out_file, "", temp_str))
            if(num == 0){
              mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
              mis_temp = mis_temp + 1
            }
          }
          
        }
      }
    }
  }
}

mis_vec = rbindlist(mis_vec_list)

#check whether the polypred converged
mis_vec_list = list()

mis_temp =  1

for(i in 1:1){
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
  
  for(l in 1:3){
    for(m in 4:4){
      for(i_rep in 1:10){
        for(i1 in 1:5){
          files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),full.names = T)
          out_file = paste0(poly_fun_file_out = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
          if(out_file%in%files==F){
            mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
            mis_temp = mis_temp + 1
          }
          
          
        }
      }
    }
  }
}

mis_vec = rbindlist(mis_vec_list)
