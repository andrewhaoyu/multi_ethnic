#goal merge gwas summary level statistics
#merge the subfiles together
eth <- c("EUR","AFR","AMR","EAS")
code <- c("c",4)
for(i in 1:4){
  temp.code <- paste0("head -1 /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary_chr1.out.assoc.linear > /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out; awk 'FNR>1'")
  for(j in 1:22){
    temp.code <- paste0(temp.code, " /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary_chr",j,".out.assoc.linear")
  }
  temp.code <- paste0(temp.code," >> /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")
  code[i] <- temp.code
}
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/merge_gwas_summary.sh",row.names = F,col.names = F,quote=F)





# for
#  awk 'FNR>1' file1 file2 file3 > bigfile
# #load all of SNPs information
# load("/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
