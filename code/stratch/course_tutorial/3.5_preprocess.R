i = 1
l = 3
i1 = 1
n.rep = 10
library(data.table)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
setwd("/data/zhangh24/multi_ethnic/")
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))

y <- y[,2+(1:n.rep),drop=F]
n <- nrow(y)
y_test_mat <- y[(100000+1):nrow(y),,drop=F]
y_out <- data.frame(ID = c(100001:120000),y = y_test_mat$V3)
write.table(y_out, file = 
              paste0("/data/BB_Bioinformatics/stat_gene_course/data/y_out"),
            row.names = F,
            col.names = T,
            quote = F)
i_rep = 1
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
summary.eur.select <- summary.eur %>% 
  filter(CHR==22)
write.table(summary.eur.select,
            file = 
              paste0("/data/BB_Bioinformatics/stat_gene_course/data/EUR_sum_data"),
            row.names = F,
            col.names = T,
            quote = F)

#to calculate PRS, we need to use PLINK --score command
#The score command needs three information: SNP, effect allele and weight
#set up the p-value threshold


subject_file <- data.frame(ID=c(100001:120000),ID2 = c(100001:120000))
write.table(subject_file, file = paste0("/data/BB_Bioinformatics/stat_gene_course/result/subject_file"),
            row.names = F,
            col.names = F,
            quote = F)

system(paste0("/data/BB_Bioinformatics/stat_gene_course/software/plink ",
              "--keep /data/BB_Bioinformatics/stat_gene_course/result/subject_file ",
              "--bfile /data/BB_Bioinformatics/stat_gene_course/data/prs_genotype/chr22 ",
              "--make-bed ",
              "--out /data/BB_Bioinformatics/stat_gene_course/data/prs_genotype/chr22_test"))