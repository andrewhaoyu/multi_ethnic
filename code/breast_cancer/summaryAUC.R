get.correlation.adj <- function(prs.model, tau, KG.plink.pre = '/data/songl5/1000G/KG.all.chr', soFile = '/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/scripts/getAdjCorrelation.so', pos_thr = 5e8){
  kg.fam = read.table(file = paste0(KG.plink.pre, '.fam'), header = F, stringsAsFactors = F)
  NSample = nrow(kg.fam)
  
  kg.bim = read.table(file = paste0(KG.plink.pre, '.bim'), header = F, stringsAsFactors = F)
  # NSNP = nrow(kg.bim)
  
  snps = intersect(prs.model[, 1], kg.bim[, 2])
  NSNP = length(snps)
  idx = match(snps, prs.model[, 1])
  prs.model = prs.model[idx, ]
  tau = tau[idx]
  
  idx = match(prs.model[, 1], kg.bim[, 2])
  chrs = kg.bim[idx, 1]
  pos = kg.bim[idx, 4]
  flip.idx = which(prs.model[, 2] != kg.bim[idx, 5])
  idx[flip.idx] = -idx[flip.idx]
  
  beta_tau = prs.model[, 3] * tau
  kg.bed.file = paste0(KG.plink.pre, '.bed')
  
  if(!file.exists(soFile)){
    cFile = gsub('.so$', '.c', soFile)
    system(paste('R CMD SHLIB', cFile))
  }
  dyn.load(soFile)
  
  adj.cor.results = .C("getAdjCorrelation", plinkBed = kg.bed.file, NumSample = as.integer(NSample), NumSNP = as.integer(NSNP), idx = as.integer(idx), chrs = as.integer(chrs), pos = as.integer(pos), pos_thr = as.integer(pos_thr), beta_tau = as.double(beta_tau), adj_cor = as.double(0.1))
  
  return(adj.cor.results$adj_cor)
}

auc <- function(prs.model.file = "/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/pruning_0.01/pv_5e-8/prs.model", gwas.summary.stats.file = "/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/MGS_0.5_1.0_CHR.SNP.A1.MAF.BETA.P.txt.gz", N0 = 2653, N1 = 2681, soFile = '/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/scripts/getAdjCorrelation.so', KG.plink.pre = '/data/songl5/1000G/KG.all.chr', pos_thr = 5e8, flag.correlation.adj.imputated.data = FALSE){
  #browser()
  library("mvtnorm")
  work.dir = dirname(prs.model.file)
  setwd(work.dir)
  
  prs.model = read.delim(file = prs.model.file, header = F, sep = '\t', stringsAsFactors = F)
  
  gwas.summary.stats = read.table(file = gwas.summary.stats.file, header = T, stringsAsFactors = F, check.names = F)
  # CHR	SNP	A1	MAF	BETA	P	INFO
  # 14	rs7160549	C	0.472562	0.206895	2.36804e-07	0.997
  snps = intersect(prs.model[, 1], gwas.summary.stats[, 'SNP'])
  
  idx = match(snps, gwas.summary.stats[, 'SNP'])
  gwas.summary.stats = gwas.summary.stats[idx, ]
  idx = match(snps, prs.model[, 1])
  prs.model = prs.model[idx, ]
  
  if("OR" %in% colnames(gwas.summary.stats))
    beta_gwas = log(gwas.summary.stats[, 'OR'])
  else
    beta_gwas = gwas.summary.stats[, 'BETA']
  idx = gwas.summary.stats[, 'A1'] != prs.model[, 2]
  beta_gwas[idx] = -beta_gwas[idx]
  P = gwas.summary.stats[, 'P']
  MAF = gwas.summary.stats[, 'MAF']
  Z = sign(beta_gwas) * qnorm(1 - P / 2)
  
  if("INFO" %in% colnames(gwas.summary.stats))
    INFO = as.numeric(gwas.summary.stats[, 'INFO'])
  else
    INFO = rep(1, length(MAF))
  
  tau = sqrt(2 * MAF * (1 - MAF) * INFO)
  beta_prs = prs.model[, 3]
  
  
  if(flag.correlation.adj.imputated.data){
    
    if(!file.exists(soFile)){
      cFile = gsub('.so$', '.c', soFile)
      system(paste('R CMD SHLIB', cFile))
    }
    dyn.load(soFile)
    
    adj.cor = 0
    for(chr in 1:22){
      cur.chr.dosage.file = paste0("dosage.chr.", chr, ".txt.gz")
      if(!file.exists(cur.chr.dosage.file))
        next
      
      # n.cur.chr.snp = as.integer(system(paste0("zcat dosage.chr.", chr, ".txt.gz | wc -l"), intern = TRUE))
      # if(n.cur.chr.snp < 1)
      # 	next
      
      cat("Chromosoe ", chr, ' Reading dosage data:\n', sep = "")
      dosage = read.delim(file = paste0("dosage.chr.", chr, ".txt.gz"), header = F, stringsAsFactors = F, sep = '\t')
      # 1     rs61769350      693731  A       G       1.944 ......
      # https://www.cog-genomics.org/plink/1.9/assoc#dosage
      # 'format=1' normally indicates a single 0..2 A1 expected count, here A1 is the 4th column with skip0=1 skip1=1 format=1
      
      cur.snps = intersect(snps, dosage[, 2])
      if(length(cur.snps) < 1)
        next
      
      idx = match(cur.snps, dosage[, 2])
      dosage = dosage[idx, ]
      dosage.chr = dosage[, 1]
      dosage.pos = dosage[, 3]
      dosage.datamatrix = t(as.matrix(dosage[, -(1:5), drop = F]))
      storage.mode(dosage.datamatrix) = 'double'
      
      idx = dosage[, 4] != prs.model[match(cur.snps, snps), 2]
      dosage.datamatrix[, idx] = 2 - dosage.datamatrix[, idx]
      
      cur.tau = rep(NA, ncol(dosage.datamatrix))
      for(i in 1:ncol(dosage.datamatrix))
        cur.tau[i] = sqrt(var(dosage.datamatrix[, i], na.rm = TRUE))
      
      idx = match(cur.snps, snps)
      cur.beta_prs = beta_prs[idx]
      cur.beta_tau = cur.beta_prs * cur.tau
      tau[idx] = cur.tau
      
      if(ncol(dosage.datamatrix) < 2){
        rm(dosage, dosage.datamatrix, dosage.chr, dosage.pos, cur.tau, cur.beta_tau)
        next
      }
      
      idx = is.na(dosage.datamatrix)
      dosage.datamatrix[idx] = -9
      
      adj.cor.results = .C("getAdjCorrelationDosage", D = as.double(dosage.datamatrix), NumSample = as.integer(nrow(dosage.datamatrix)), NumSNP = as.integer(ncol(dosage.datamatrix)), chrs = as.integer(dosage.chr), pos = as.integer(dosage.pos), pos_thr = as.integer(pos_thr), beta_tau = as.double(cur.beta_tau), adj_cor = as.double(0.1))
      adj.cor = adj.cor + adj.cor.results$adj_cor
      
      rm(dosage, dosage.datamatrix, dosage.chr, dosage.pos, cur.tau, cur.beta_tau, adj.cor.results)
    }
  }else{
    
    adj.cor = get.correlation.adj(prs.model, tau, KG.plink.pre, soFile, pos_thr)
    # cat(paste0('adj.cor:', adj.cor, '\n'))
    # delta = sqrt(1/N1 + 1/N0) * sum(beta_prs * Z * tau) / sqrt(2 * sum((beta_prs * tau)^2) + 4 * adj.cor)
    # cat(paste0('delta:', delta, '\n'))
  }
  
  # cat(paste0('adj.cor:', adj.cor, '\n'))
  
  delta = sqrt(1/N1 + 1/N0) * sum(beta_prs * Z * tau) / sqrt(2 * sum((beta_prs * tau)^2) + 4 * adj.cor)
  # cat(paste0('delta:', delta, '\n'))
  
  auc0 = pnorm(delta)
  pxy = pmvnorm(lower=c(-delta,-delta),upper=Inf,mean=c(0,0),sigma=matrix(c(1,0.5,0.5,1),2,2))[1]
  var0 = (pxy - auc0^2) * (N1+N0) / (N1*N0)
  return(c(auc0, var0))
}


res = auc(	prs.model.file = "/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/pruning_0.01/pv_5e-8/prs.model", 
           gwas.summary.stats.file = "/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/MGS_0.5_1.0_CHR.SNP.A1.MAF.BETA.P.txt.gz",
           N0 = 2653,
           N1 = 2681,
           soFile = '/data/songl5/Jianxin/AUC.SCZ/data/pruning_imputation_withPCA/pruning_0.5_1.0/scripts/getAdjCorrelation.so',
           flag.correlation.adj.imputated.data = FALSE,
           pos_thr = 5e8,
           KG.plink.pre = '/data/songl5/1000G/KG.all.chr')

cat("\n#######################################\n\n")

cat("Predicted AUC:\t", res[1], "\n", sep = "")
cat("Predicted AUC's variance:\t", res[2], "\n", sep = "")

cat("\nHave a nice Day!\n")

cat("\n#######################################\n")
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
