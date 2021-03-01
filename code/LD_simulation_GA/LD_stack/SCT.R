#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
i_rep = as.numeric(args[[3]])
#l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
#i_rep = 1

#j = as.numeric(args[[3]])
library(dplyr)
library(data.table)
library(bigsnpr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n.snp.mat <- matrix(0,length(pthres),4)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir,"."))
system(paste0("ls ",temp.dir))
bim <- as.data.frame(fread(paste0(temp.dir,"all_chr_test.mega.bim")))
colnames(bim)[2] <- "SNP"
print("step1 finished")
snp_readBed(paste0(temp.dir,"all_chr_test.mega.bed"))
obj.bigSNP <- snp_attach(paste0(temp.dir,"all_chr_test.mega.rds"))
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# See how the file looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
big_counts(G, ind.col = 1:10)
NCORES <- 2
#for(l in 1:3){

setwd("/data/zhangh24/multi_ethnic/")
library(tidyverse)
#for(m in 1:4){
  sumstats <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
  sumstats = left_join(bim,sumstats,by="SNP") %>% 
    distinct(SNP,.keep_all = T)
  
  #snp.id = sumstats %>%  select(SNP)
  #get the reference allele and effect allele
  #sum data used A1 to code effect allele
  #bim data used original plink coding
  #bigsnpr required both reference and effect
  #transform sum data to orginal plink coding
  eff_allele = sumstats$A1
  noneff_allele = sumstats$V6
  coding_allele = sumstats$V5
  noncoding_allele = sumstats$V6
  BETA = sumstats$BETA
  idx <- which(eff_allele!=coding_allele)
  eff_allele[idx] = coding_allele[idx]
  BETA[idx] = -BETA[idx]
  noneff_allele[idx] = noncoding_allele[idx]
  all.equal(eff_allele,coding_allele)
  sumstats$a0 = eff_allele
  sumstats$a1 = noneff_allele
  sumstats$BETA = BETA
  
  sumstats = sumstats %>% 
    rename(chr=CHR,
           rsid = SNP,
           pos = BP,
           beta = BETA,
           p = P) %>% 
    select(chr,rsid,pos,a0,a1,beta,p)
  # ind.train <- c(1:10000)
  # ind.test <- setdiff(rows_along(G), ind.train)
  # 
  ind.train <-  c(1:3000)
  ind.test <- setdiff(rows_along(G), ind.train)
  map <- obj.bigSNP$map[,-(2:3)]
  names(map) <- c("chr", "pos", "a0", "a1")
  info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
  info.snp.id = info_snp$rsid
  info.snp.chr = info_snp$chr
  info.snp.pos = info_snp$pos
  beta <- info_snp$beta
  lpval <- -log10(info_snp$p)
  all_keep <- snp_grid_clumping(G, info.snp.chr, info.snp.pos, ind.row = ind.train,
                                lpS = lpval, ncores = NCORES)
  ind.train = c(1:10000)
  multi_PRS = snp_grid_PRS(G, all_keep, beta, lpval, ind.row = ind.train,
                           backingfile = paste0(temp.dir,"all_chr_test.mega.rds"), 
                           n_thr_lpS = 50, ncores = NCORES)
  phe.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  y <- as.data.frame(fread(paste0(phe.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
  n.rep = 10
  y <- y[,2+(1:n.rep),drop=F]
  n <- nrow(y)
  y_test_mat <- y[100001:nrow(y),,drop=F]
  y = y_test_mat[,i_rep]
  final_mod <- snp_grid_stacking(multi_PRS, y[ind.train], ncores = NCORES, K = 10)
  summary(final_mod$mod)
  new_beta <- final_mod$beta.G
  ind <- which(new_beta != 0)
  ind.test = c(10001:20000)
  pred <- 
    big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)
  y.vad = y_test_mat[10001:20000,i_rep]
  model2 <- lm(y.vad~pred)
  r2 = summary(model2)$r.square
  save(r2,file = paste0(out.dir,eth[i],"/r2.SCT_rho_",l,"_size_",m,"_GA_",i1,"_rep_",i_rep))
 




# result.matrix <- matrix(0,3,1)
# for(l in 1:3){
#   #for(i1 in 1:2){
#     sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep)))
#     idx <- which(sum.data$P<=5E-08)
# 
#     result.matrix[l,i1] <- length(idx)
#   #}
# }
