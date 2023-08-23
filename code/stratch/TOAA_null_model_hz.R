library(devtools)
devtools::install_github("xihaoli/STAAR")
library(STAAR)
TOAA_phenotype <- read.csv("/data/COVID_ADHOC/lxwg/STAARpipeline/Steps/Association_Analysis/Step1/Input/TOAA_Phenotype_10_PCs.csv")
# fit null model assuming subjects are unrelated
obj.STAAR.unrelated <- glm(trait~age+as.factor(gender)+scale(PC1)+
                             scale(PC2)+scale(PC3)+scale(PC4)+scale(PC5), data = TOAA_phenotype, family = binomial(link = "logit"))
# Warning messages:
# 1: glm.fit: fitted probabilities numerically 0 or 1 occurred

sgrm <- get(load("/data/COVID_ADHOC/lxwg/STAARpipeline/Steps/Association_Analysis/Step1/Input/output.sparseGRM.sGRM.RData"))
sample_id <- unlist(lapply(colnames(sgrm),function(x){substring(x,1,floor(nchar(x)/2))}))
colnames(sgrm) <- sample_id
rownames(sgrm) <- sample_id

### fit null model assuming subjects are related
obj.STAAR <- fit_null_glmmkin(trait~age+gender+scale(PC1)+
                                scale(PC2)+scale(PC3), data = TOAA_phenotype, kins = sgrm, kins_cutoff = 0.022, id = "SampleID", use_sparse = TRUE,family = binomial(link = "logit"), verbose=T)
