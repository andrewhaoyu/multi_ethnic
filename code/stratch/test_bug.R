#use GCTA to generate phenotypes
#generate phenotypes given different genetic architecture
args = commandArgs(trailingOnly = T)
#i represent ethnic group
i = 1
#l represent causal snps proportion
l = 3
#i1 represent genetic architecture
i1 = 2
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"


# sid<-Sys.getenv('SLURM_JOB_ID')
# dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
eth <- c("EUR","AFR","AMR","EAS","SAS")
# system(paste0("cp ",cur.dir,eth[i],"/select.cau.snp.bed /lscratch/",sid,"/test/",eth[i],"_select.cau.snp.bed"))
# system(paste0("cp ",cur.dir,eth[i],"/select.cau.snp.bim /lscratch/",sid,"/test/",eth[i],"_select.cau.snp.bim"))
# system(paste0("cp ",cur.dir,eth[i],"/select.cau.snp.fam /lscratch/",sid,"/test/",eth[i],"_select.cau.snp.fam"))

select.cau <- read.table(paste0(out.dir,eth[i],"/select.cau_rho",l,"_",i1),header=F)

snp.infor <- read.table("/data/zhangh24/test_code/select.cau.snp.bim",header=F)
idx <- which(snp.infor$V2=="rs193160839:740285:G:A")
snp.infor[idx,]

snp.freq <- read.table("/data/zhangh24/test_code/freq_1.frq",header=T)

res <- system(paste0("/data/zhangh24/software/plink2 --bfile  /data/zhangh24/test_code/select.cau.snp --freq --out  /data/zhangh24/test_code/freq_1"))

head(select.cau)
jdx <- which(select.cau$V1=="rs12137148:2768649:A:G")
select.cau[jdx,]

select.cau <- select.cau[1,]
select.cau[,2] <- -0.05
#herit <- nrow(select.cau)*var(select.cau$V2)
write.table(select.cau,file = "/data/zhangh24/test_code/select.cau.txt",row.names = F,col.names = F,quote = F)

res <- system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile  /data/zhangh24/test_code/select.cau.snp --simu-qt --simu-causal-loci /data/zhangh24/test_code/select.cau.txt --simu-hsq 1 --simu-rep 1 --out  /data/zhangh24/test_code/phenotype_1"))
pheno1 <- read.table("/gpfs/gsfs11/users/zhangh24/test_code/phenotype_1.phen")
head(pheno1)
table(pheno1$V3)
out.put <- read.table("/data/zhangh24/test_code/phenotype_1.par",header=T)
head(out.put)
range(out.put$Frequency)
jdx <- which(out.put$QTL=="rs12137148:2768649:A:G")
out.put[jdx,]








snp.infor.alter <- snp.infor
snp.infor.alter[idx,6] <- snp.infor[idx,5]
snp.infor.alter[idx,5] <- snp.infor[idx,6]
snp.infor.alter[idx,]
write.table(snp.infor.alter,file = "/data/zhangh24/test_code/select.cau.snp.alter.bim",row.names = F,col.names = F,quote=F)
res <- system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile  /data/zhangh24/test_code/select.cau.snp.alter --simu-qt --simu-causal-loci /data/zhangh24/test_code/select.cau.txt --simu-hsq 1 --simu-rep 1 --out  /data/zhangh24/test_code/phenotype_2"))
pheno2 <- read.table("/gpfs/gsfs11/users/zhangh24/test_code/phenotype_2.phen")
head(pheno2)
table(pheno2$V3)

select.cau.test <- select.cau
RefAllele <- "C"
Frequency <- 0.017
select.cau.test <- data.frame(select.cau.test,RefAllele,Frequency)
select.cau.test <- select.cau.test[,c(1,3,4,2)]
colnames(select.cau.test) <- c("QTL","RefAllele","Frequency","Effect")
write.table(select.cau.test,file = "/data/zhangh24/test_code/select.cau.test.txt",row.names = F,col.names = T,quote = F)
res <- system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile  /data/zhangh24/test_code/select.cau.snp.alter --simu-qt --simu-causal-loci /data/zhangh24/test_code/select.cau.test.txt --simu-hsq 1 --simu-rep 1 --out  /data/zhangh24/test_code/phenotype_3"))
pheno3 <- read.table("/gpfs/gsfs11/users/zhangh24/test_code/phenotype_3.phen")
head(pheno3)
table(pheno3$V3)

