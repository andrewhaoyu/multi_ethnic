library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]
newpsa$l_marker1 <- log(newpsa$marker1)
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
                     group = "status",
                     tag.h = 0,
                     data = newpsa,
                     standardise = TRUE,
                     p = seq(0,1,l=101),
                     compute.lpml = TRUE,
                     compute.WAIC = TRUE,
                     compute.DIC = TRUE,
                     pauc = pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
                     density = densitycontrol.aroc(compute = TRUE, grid.h = NA, newdata = NA),
                     prior.h = priorcontrol.bnp(m0 = rep(0, 4), S0 = 10*diag(4), nu = 6, Psi = diag(4),
                                                a = 2, b = 0.5, alpha = 1, L =10),
                     mcmc = mcmccontrol(nsave = 500, nburn = 100, nskip = 1))
tmp = summary(AROC_bnp)
as.numeric(strsplit(gsub("Area under the covariate-adjusted ROC curve:","",tmp$AUC)," ")[[1]][[2]])
  



auc.result.vec = rep(0,100)
auc.result.sl.vec = rep(0,100)
for(si in 1:100){
  load(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/auc_tdld_eb",si,".rdata"))  
  auc.result.vec[si] = auc.tdld.eb[[1]]
  auc.result.sl.vec[si] = auc.tdld.eb[[3]]
}
which.max(auc.result.vec)
which.max(auc.result.sl.vec)

idx <- which.max(result.data$auc.test.vec)

library(data.table)
data1 = fread("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/EUR/summary_combine/summary_rho_1_size_1_GA_1")
data2 = fread("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/EUR/summary_combine/summary_rho_1_size_1_GA_5")