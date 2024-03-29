#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
i_rep = as.numeric(args[[1]])
i = as.numeric(args[[2]])
#j = as.numeric(args[[3]])
l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
#i_c=1 for mega; i_c=2 for 1kg
i1 = 1
i_c = as.numeric(args[[5]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/test_chip/"
#j = as.numeric(args[[3]])
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
if(i_c==2){
  #hapmap snps just need the mega subset
  system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))
  system(paste0("ls ",temp.dir))
  
}else if(i_c==1){
  #1kg snps
  
  system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.tag.* ",temp.dir))
}

library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/")

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
#n.snp.mat <- matrix(0,length(pthres),4)


#for(l in 1:3){
print(l)
# for(m in 1:4){

print(m)
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>% 
  select(SNP,beta_eur,peur)
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com <- left_join(sum.data,summary.eur.select,by="SNP")
#remove duplicated snp id
idx <- which(duplicated(summary.com$SNP))
dup.id = summary.com$SNP[idx]
summary.com = summary.com %>% filter(SNP%in%dup.id==F)

rs.id.list = list()
temp = 1
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_i_c_",i_c,"rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    rs.id.list[[temp]] = LD   
    temp = temp + 1
  }
}
unique.id = unique(rbindlist(rs.id.list))
write.table(unique.id,file = paste0(temp.dir,"unique.id"),row.names = F,col.names = T,quote=F)
#extract the subset of SNPs used in prs calculation
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--bfile ",temp.dir,"all_chr_test.tag ",
                    "--extract ",temp.dir,"unique.id ",
                    "--make-bed ",
                    "--out ",temp.dir,"extact_test"))  
system(paste0("rm ",temp.dir,"all_chr_test.tag*"))








for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    
    
    #read LD clumped SNPs
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_i_c_",i_c,"rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    
    #read the target ethnic group summary level statistics
    #combine the statistics with SNPs after clumping
    prs.all <- left_join(LD,summary.com,by="SNP") 
    colSums(is.na(prs.all))
    
    for(k1 in 1:length(pthres)){
      #keep al the SNPs with peur pass the threshold
      
      
      prs.file = prs.all %>% 
        mutate(P = replace(P,peur<=pthres[k1],1E-20)) %>% 
        select(SNP,A1,BETA,P) %>% 
        mutate(BETA_2 = rnorm(nrow(prs.file),0))
      write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
      
      p.value.file <- prs.file %>% 
        select(SNP,P)
      write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
      n_pthres = length(pthres)
      q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
      temp = 1
      
      for(k2 in 1:length(pthres)){
        prs.file.sub <- prs.file %>% 
          filter(P<=pthres[k2]) %>% 
          select(SNP,A1,BETA,P)
        if(nrow(prs.file.sub)>0){
          q_range[temp,1] = paste0("p_value_",k2)
          q_range[temp,3] = pthres[k2]
          temp = temp+1
        }
      }
      q_range = q_range[1:(temp-1),]
      write.table(q_range,file = paste0(temp.dir.prs,"q_range_file"),row.names = F,col.names = F,quote=F)
      if(i_c==2){
        res = system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"all_chr_test.mega --exclude ",old.out.dir,eth[i],"/duplicated.id  --out ",temp.dir.prs,"prs_2DLD_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))  
        
        
      }else{
        # fam <- as.data.frame(fread(paste0(temp.dir,eth[i],"/all_chr.tag.fam")))
        # ref_fam = fam[100001:120000,]
        # write.table(ref_fam,file = paste0(temp.dir,"/ref_fam.fam"),row.names=F,col.names=F,quote=F)
        # 
    
        res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                            "--q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header ",
                            "--threads 2 ",
                            "--score-col-nums 3,5 ",
                            "--score ",temp.dir.prs,"prs_file header no-mean-imputation ",
                            "--bfile ",temp.dir,"extact_test ",
                            "--exclude ",old.out.dir,eth[i],"/duplicated.id  ",
                            "--out ",temp.dir.prs,"prs_2DLD_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))  
    system(paste0("ls /lscratch/32148155/test/prs/prs_2DLD_rho_1_size_1_rep_1_GA_1_rind_1_wcind_1p_value_8.*.sscore"))
    
    
        }
      
      print("step2 finished")
      #system(paste0("ls ",temp.dir.prs))
      if(res==2){
        stop()
      }
      
      
      
      #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
    }
    
  }
  
}

total = length(r2_vec)*length(wc_vec)*length(pthres)*length(pthres)
prs.list = list()
name.file = rep("c",total)

temp = 1
#files = dir(temp.dir.prs,pattern="prs_2DLD_rho_",full)
if(i_c==2){
  #hapmap used plink 1.9
  for(r_ind in 1:length(r2_vec)){
    wc_vec = wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      for(k1 in 1:length(pthres)){
        for(k2 in 1:length(pthres)){
          filename = paste0(temp.dir.prs,"prs_2DLD_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
          #if(filename%in%files){
          temp.prs.file = fread(filename)
          prs.list[[temp]] = temp.prs.file[,6,drop=F]
          name.file[temp] = paste0("r_ind_",r_ind,"_w_ind_",w_ind,"_k1_",k1,"_k2_",k2)
          temp = temp + 1
          
          # }
        }
      }
    }
    
  }
  name.file = name.file[1:(temp-1)]
  prs.mat = bind_cols(prs.list)
  colnames(prs.mat) = name.file
  prs.id = temp.prs.file[,1:2,drop=F]
  prs.file = cbind(prs.id,prs.mat)
  write.table(prs.file,file = paste0(out.dir,eth[i],"/prs/prs_2DLD_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".profile"))
  
}else{
  #1kg used plink 2.0 for speed
  for(r_ind in 1:length(r2_vec)){
    wc_vec = wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      for(k1 in 1:length(pthres)){
        for(k2 in 1:length(pthres)){
          filename = paste0(temp.dir.prs,"prs_2DLD_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".sscore")
          #if(filename%in%files){
          temp.prs.file = fread(filename)
          prs.list[[temp]] = temp.prs.file[,5,drop=F]
          name.file[temp] = paste0("r_ind_",r_ind,"_w_ind_",w_ind,"_k1_",k1,"_k2_",k2)
          temp = temp + 1
          
          # }
        }
      }
    }
    
  }
  name.file = name.file[1:(temp-1)]
  prs.mat = bind_cols(prs.list)
  colnames(prs.mat) = name.file
  prs.id = temp.prs.file[,1:2,drop=F]
  prs.file = cbind(prs.id,prs.mat)
  write.table(prs.file,file = paste0(out.dir,eth[i],"/prs/prs_2DLD_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".profile"))
  
}



#data = fread(paste0(out.dir,eth[i],"/prs/prs_2DLD_i_c_",i_c,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".profile"))
