#goal merge gwas summary level statistics
#merge the subfiles together


cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"

#five ethnics i
eth <- c("EUR","AFR","AMR","EAS","SAS")
#three different causal proportion l
#three different training sample sizes m
#22 chr
code <- c("c",5*3*3*10)
temp = 1
for(i in 1:5){
  for(l in 1:3){
    for(i_rep in 1:10){
      for(m in 1:4){
        temp.code <- paste0("head -1 ",cur.dir,eth[i],"/summary_chr_",1,"_rho_",l,"_rep_",i_rep,".out.P",m,".assoc.linear > ",cur.dir,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"; awk 'FNR>1'")
        for(j in 1:22){
          temp.code <- paste0(temp.code, " ",cur.dir,eth[i],"/summary_chr_",j,"_rho_",l,"_rep_",i_rep,".out.P",m,".assoc.linear")
        }
        temp.code <- paste0(temp.code," >> ",cur.dir,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep)
        code[temp] <- temp.code
        temp <- temp+1
      }
    }
  }  
    }
    
 
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation_large/merge_gwas_summary.sh",row.names = F,col.names = F,quote=F)

      



# for
#  awk 'FNR>1' file1 file2 file3 > bigfile
# #load all of SNPs information
# load("/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
