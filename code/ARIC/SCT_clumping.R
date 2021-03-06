#LD_clumping for ARIC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

#i = as.numeric(args[[1]])
#l = as.numeric(args[[2]])
j = as.numeric(args[[1]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
library(bigsnpr)
for(i in 1:2){
  for(l in c(1,3)){
    setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    #load gwas summary statistics
    sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    # write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
    #             ,col.names = T,row.names = F,quote=F) 
    #prepare association file for plink
    sum.data.assoc = sum.data %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL) %>% 
      filter(CHR==j) 
    
    sum.data.assoc = sum.data.assoc[,c("CHR","SNP","BP","A1","BETA","P")]
    #snp_readBed(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/ref_chr",j,".bed"))
    obj.bigSNP <- snp_attach(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j,".rds"))
    str(obj.bigSNP, max.level = 2, strict.width = "cut")
    # See how the file looks like
    str(obj.bigSNP, max.level = 2, strict.width = "cut")
    G   <- obj.bigSNP$genotypes
    CHR <- obj.bigSNP$map$chromosome
    POS <- obj.bigSNP$map$physical.pos
    y   <- obj.bigSNP$fam$affection - 1
    big_counts(G, ind.col = 1:20)
    NCORES <- 2
    map <- obj.bigSNP$map
    colnames(map)[2] = "SNP"
    #get the reference allele and effect allele
    #sum data used A1 to code effect allele
    #bim data used original plink coding
    #map file minor allele is A1 major allele is A2. T
    #The Genotype data in bed file is coded with minor allele A1
    #bigsnpr required both reference and effect
    #transform sum data to orginal plink coding
    sum.data.match <- left_join(map,sum.data.assoc)
    
    #idx <- which(sum.data.match$A1!=sum.data.match$allele1)
    sum.data.match = sum.data.match %>% 
      mutate(effect_allele = A1,
             noneffect_allele = allele2) %>% 
      mutate(BETA = ifelse(A1==allele1,BETA,-BETA),
             a0 = allele1,
             a1 = allele2) %>% 
      rename(chr = CHR,
             pos = physical.pos,
             beta = BETA) %>% 
      select(chr,SNP,pos,a0,a1,beta,P)
    
    map <- obj.bigSNP$map[,-(2:3)]
    names(map) <- c("chr", "pos", "a0", "a1")
    
    
    
    info_snp <- snp_match(sum.data.match, map, strand_flip = FALSE)
    beta <- rep(NA, ncol(G))
    beta[info_snp$`_NUM_ID_`] <- info_snp$beta
    lpval <- rep(NA, ncol(G))
    lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)
    ind.train = sample(nrow(G), 1000)
    all_keep <- snp_grid_clumping(G, CHR, POS,ind.row = ind.train,
                                  lpS = lpval, exclude = which(is.na(lpval)),ncores = NCORES)
    save(all_keep,file = paste0(temp.dir,"all_keep_chr_",j,".rdata"))
  }
}
    #idx <- which(sum.data.assoc$SNP=="rs4970836")
    #write.table(sum.data.assoc,file = paste0(temp.dir,"chr_",j,"_assoc.txt"),col.names = T,row.names = F,quote=F)
    
    # dim(summary)
    # head(summary)
  #   pthr = 0.5
  #   r2thr = 0.1
  #   kbpthr = 500
  #   eth <- c("EUR","AFR","AMR","EAS","SAS")
  #   #cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  #   #code <- rep("c",5*3*3)
  #   #system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
  #   res = system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/ref_chr",j," --clump ",temp.dir,"chr_",j,"_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_clump_chr_",j))
  #   system(paste0("mv ",temp.dir,"LD_clump_chr_",j,".clumped ",out.dir))
  #   if(res==2){
  #     stop()
  #   }
  #   
  #   paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/ref_chr",j)
  # }
}
