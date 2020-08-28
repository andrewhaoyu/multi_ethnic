#goal remove extra files
args = commandArgs(trailingOnly = T)

i = as.numeric(args[[1]])
#j = as.numeric(args[[2]])
#l = as.numeric(args[[3]])
#m = as.numeric(args[[4]])



cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"

#five ethnics i
eth <- c("EUR","AFR","AMR","EAS","SAS")
#three different causal proportion l


        for(i_rep in 11:100){
          system(paste0("rm ",cur.dir,eth[i],"/summary_chr_*_rho_*_rep_",i_rep,".out.P*.assoc.linear"))  
        }
        
