#Goal: select the causal SNPs
#the heritability are attributable to different proportion of SNPs
#the total heritability are the same
#different proportion of causal SNPs

#load the SNPs information
load("/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
#only select SNPs that don't have duplicated position number
#It will make the extract procedure easier
#the CHR2 EUR file after 700000 are corrupted
#causal SNPs shouldn't come from there
# idx <- which(snp.infor$EUR>=0.01&
#                snp.infor$EUR<=0.99&
#                snp.infor$CHR==2&
#                snp.infor$position>=239554914)
# snp.infor <- snp.infor[-idx,]
# save(snp.infor,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")

EUR <- which(snp.infor$EUR>=0.01&
               snp.infor$EUR<=0.99)
AFR <-  which(snp.infor$AFR>=0.01&
                snp.infor$AFR<=0.99)
AMR <-  which(snp.infor$AMR>=0.01&
                snp.infor$AMR<=0.99)
EAS <-  which(snp.infor$EAS>=0.01&
                snp.infor$EAS<=0.99)



n.cau.EUR <- 5000
n.cau.AFR <- round(n.cau.EUR*length(AFR)/length(EUR),0)
#due to round error, just + 1
n.cau.AMR <- round(n.cau.EUR*length(AMR)/length(EUR),0)
n.cau.EAS <- round(n.cau.EUR*length(EAS)/length(EUR),0)
set.seed(666)
all.id <- c(1:nrow(snp.infor))


#group1 EUR only
EUR.only <- which(all.id%in%EUR==T&
                    all.id%in%AFR==F&
                    all.id%in%AMR==F&
                    all.id%in%EAS==F)
n.cau.EUR.only <- round(n.cau.EUR*length(EUR.only)/length(EUR),0)
cau.EUR.only <- sample(EUR.only,n.cau.EUR.only,replace=F)
#group2 EUR AFR 
EUR.AFR <- which(all.id%in%EUR==T&
                   all.id%in%AFR==T&
                   all.id%in%AMR==F&
                   all.id%in%EAS==F)
n.cau.EUR.AFR <- round(n.cau.EUR*length(EUR.AFR)/length(EUR),0)
cau.EUR.AFR <- sample(EUR.AFR,n.cau.EUR.AFR,replace=F)
#group3 EUR AMR
EUR.AMR <- which(all.id%in%EUR==T&
                   all.id%in%AFR==F&
                   all.id%in%AMR==T&
                   all.id%in%EAS==F)
n.cau.EUR.AMR <- round(n.cau.EUR*length(EUR.AMR)/length(EUR),0)+1
cau.EUR.AMR <- sample(EUR.AMR,n.cau.EUR.AMR,replace=F)
#group4 EUR EAS
EUR.EAS <- which(all.id%in%EUR==T&
                   all.id%in%AFR==F&
                   all.id%in%AMR==F&
                   all.id%in%EAS==T)
n.cau.EUR.EAS <- round(n.cau.EUR*length(EUR.EAS)/length(EUR),0)
cau.EUR.EAS <- sample(EUR.EAS,n.cau.EUR.EAS,replace=F)
#group5 EUR AFR AMR
EUR.AFR.AMR <- which(all.id%in%EUR==T&
                       all.id%in%AFR==T&
                       all.id%in%AMR==T&
                       all.id%in%EAS==F)
n.cau.EUR.AFR.AMR <- round(n.cau.EUR*length(EUR.AFR.AMR)/length(EUR),0)
cau.EUR.AFR.AMR <- sample(EUR.AFR.AMR,n.cau.EUR.AFR.AMR,replace=F)
#group6 EUR AFR EAS
EUR.AFR.EAS <- which(all.id%in%EUR==T&
                       all.id%in%AFR==T&
                       all.id%in%AMR==F&
                       all.id%in%EAS==T)
n.cau.EUR.AFR.EAS <- round(n.cau.EUR*length(EUR.AFR.EAS)/length(EUR),0)
cau.EUR.AFR.EAS <- sample(EUR.AFR.EAS,n.cau.EUR.AFR.EAS,replace=F)

#group7 EUR AMR EAS
EUR.AMR.EAS <- which(all.id%in%EUR==T&
                       all.id%in%AFR==F&
                       all.id%in%AMR==T&
                       all.id%in%EAS==T)
n.cau.EUR.AMR.EAS <- round(n.cau.EUR*length(EUR.AMR.EAS)/length(EUR),0)
cau.EUR.AMR.EAS <- sample(EUR.AMR.EAS,n.cau.EUR.AMR.EAS,replace=F)

#group8 EUR AFR AMR EAS
EUR.AFR.AMR.EAS <- which(all.id%in%EUR==T&
                           all.id%in%AFR==T&
                           all.id%in%AMR==T&
                           all.id%in%EAS==T)
n.cau.EUR.AFR.AMR.EAS <- round(n.cau.EUR*length(EUR.AFR.AMR.EAS)/length(EUR),0)
cau.EUR.AFR.AMR.EAS <- sample(EUR.AFR.AMR.EAS,n.cau.EUR.AFR.AMR.EAS,replace=F)

cau.EUR = c(cau.EUR.only,
            cau.EUR.AFR,
            cau.EUR.AMR,
            cau.EUR.EAS,
            cau.EUR.AFR.AMR,
            cau.EUR.AFR.EAS,
            cau.EUR.AMR.EAS,
            cau.EUR.AFR.AMR.EAS)
length(cau.EUR)
#group9 AFR only
AFR.only<- which(all.id%in%EUR==F&
                   all.id%in%AFR==T&
                   all.id%in%AMR==F&
                   all.id%in%EAS==F)
n.cau.AFR.only <- round(n.cau.AFR*length(AFR.only)/length(AFR),0)
cau.AFR.only <- sample(AFR.only,n.cau.AFR.only,replace=F)
#group10 AFR,AMR 
AFR.AMR <- which(all.id%in%EUR==F&
                   all.id%in%AFR==T&
                   all.id%in%AMR==T&
                   all.id%in%EAS==F)
n.cau.AFR.AMR <- round(n.cau.AFR*length(AFR.AMR)/length(AFR),0)
cau.AFR.AMR <- sample(AFR.AMR,n.cau.AFR.AMR,
                      replace=F)
#group11 AFR,EAS
AFR.EAS <- which(all.id%in%EUR==F&
                   all.id%in%AFR==T&
                   all.id%in%AMR==F&
                   all.id%in%EAS==T)
n.cau.AFR.EAS <- round(n.cau.AFR*length(AFR.EAS)/length(AFR),0)
cau.AFR.EAS <- sample(AFR.EAS,n.cau.AFR.EAS,
                      replace=F)
#group12 AFR,AMR,EAS
AFR.AMR.EAS <- which(all.id%in%EUR==F&
                       all.id%in%AFR==T&
                       all.id%in%AMR==T&
                       all.id%in%EAS==T)
n.cau.AFR.AMR.EAS <- round(n.cau.AFR*length(AFR.AMR.EAS)/length(AFR),0)
cau.AFR.AMR.EAS <- sample(AFR.AMR.EAS,n.cau.AFR.AMR.EAS,
                          replace=F)
cau.AFR = c(cau.AFR.only,
            cau.EUR.AFR,
            cau.AFR.AMR,
            cau.AFR.EAS,
            cau.EUR.AFR.AMR,
            cau.EUR.AFR.EAS,
            cau.AFR.AMR.EAS,
            cau.EUR.AFR.AMR.EAS)
length(cau.AFR)
#group13 AMR only
AMR.only <- which(all.id%in%EUR==F&
                    all.id%in%AFR==F&
                    all.id%in%AMR==T&
                    all.id%in%EAS==F)
n.cau.AMR.only <- round(n.cau.AMR*length(AMR.only)/length(AMR),0)
cau.AMR.only <- sample(AMR.only,n.cau.AMR.only,
                       replace=F)

#group14 AMR.EAS
AMR.EAS <- which(all.id%in%EUR==F&
                   all.id%in%AFR==F&
                   all.id%in%AMR==T&
                   all.id%in%EAS==T)
n.cau.AMR.EAS <- round(n.cau.AMR*length(AMR.EAS)/length(AMR),0)
cau.AMR.EAS <- sample(AMR.EAS,n.cau.AMR.EAS,
                      replace=F)

cau.AMR = c(cau.AMR.only,
            cau.EUR.AMR,
            cau.AFR.AMR,
            cau.AMR.EAS,
            cau.EUR.AFR.AMR,
            cau.EUR.AMR.EAS,
            cau.AFR.AMR.EAS,
            cau.EUR.AFR.AMR.EAS)
length(cau.AMR)
n.cau.AMR
#group15 EAS only
EAS.only <- which(all.id%in%EUR==F&
                    all.id%in%AFR==F&
                    all.id%in%AMR==F&
                    all.id%in%EAS==T)
n.cau.EAS.only <- round(n.cau.EAS*length(EAS.only)/length(EAS),0)
cau.EAS.only <- sample(EAS.only,n.cau.EAS.only,
                       replace=F)

cau.EAS= c(cau.EAS.only,
           cau.EUR.EAS,
           cau.AFR.EAS,
           cau.AMR.EAS,
           cau.EUR.AFR.EAS,
           cau.EUR.AMR.EAS,
           cau.AFR.AMR.EAS,
           cau.EUR.AFR.AMR.EAS)
length(cau.EAS)
n.cau.EAS

cau.SNP <- unique(c(cau.EUR,cau.AFR,cau.AMR,cau.EAS))

cau.snp.infor <- snp.infor[cau.SNP,]
save(cau.snp.infor,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/cau.snp.infor.rdata")
