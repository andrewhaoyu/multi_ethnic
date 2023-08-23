library(devtools)
devtools::install_github("UW-GAC/GENESIS")
devtools::install_github("xihaoli/STAAR")
library(STAAR)
TOAA_phenotype <- read.csv("~/Desktop/temp4/input/TOAA_Phenotype.csv")
#test the correlation between pc and trait
result = rep(0,10)
for(i in 1:10){
  cor.result = cor.test(TOAA_phenotype$trait,scale(TOAA_phenotype[,i+4]))
  result[i] = cor.result$p.value
}
#find which PCs are correlated with the outcome
which(result<=0.05/10)
#turn out to be pc1-3,pc5,pc6,pc7,pc8

# fit null model assuming subjects are unrelated
obj.STAAR.unrelated <- glm(as.factor(trait)~age+as.factor(gender)+scale(PC1)+
                             +scale(PC2)+scale(PC3)+scale(PC4)
                             +scale(PC5), data = TOAA_phenotype, family = binomial(link = "logit"))
summary(obj.STAAR.unrelated)

obj.STAAR.unrelated <- glm(as.factor(trait)~age+as.factor(gender)+scale(PC1)+
                             +scale(PC2)+scale(PC3)+scale(PC4)
                           +scale(PC5)+scale(PC6), data = TOAA_phenotype, family = binomial(link = "logit"))
summary(obj.STAAR.unrelated)


table(TOAA_phenotype$trait,TOAA_phenotype$gender)
table(TOAA_phenotype$trait,TOAA_phenotype$age)
summary(obj.STAAR.unrelated)
# Warning messages:
# 1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred

sgrm <- get(load("~/Desktop/Input/output.sparseGRM.sGRM.RData"))
sample_id <- unlist(lapply(colnames(sgrm),function(x){substring(x,1,floor(nchar(x)/2))}))
colnames(sgrm) <- sample_id
rownames(sgrm) <- sample_id

### fit null model assuming subjects are related
obj.STAAR <- fit_null_glmmkin(trait~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = TOAA_phenotype, kins = sgrm, kins_cutoff = 0.022, id = "SampleID", use_sparse = TRUE,family = binomial(link = "logit"), verbose=T)
var(scale(TOAA_phenotype$PC1))


library(ggplot2)
ggplot(data = TOAA_phenotype) +
  geom_point(aes(x = PC1,y = PC2, col =as.factor(trait )))

             