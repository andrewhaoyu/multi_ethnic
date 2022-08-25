# compute_r2.R computes and prints R2.
#
# The variable file_path should match with file_out in run_XPASS.R.

## clear workspace
rm(list = ls())

#####
# Install packages
packages <- c("data.table")
to_be_installed <- setdiff(packages, rownames(installed.packages()))
if (length(to_be_installed)) {
  install.packages(to_be_installed)
}

#####
# Load library
library(data.table)

#####
# input parameters
race = "AFR" # target population
size = 1 # size for target population
GA = 1
rho = 1
rep = 1
file_path = paste0("/dcs04/nilanjan/data/wlu/XPASS/XPASS_result/", race, "/XPASS-rho",rho,'-size',size,'-rep',rep,'-GA',GA)
file_path = "/users/wlu/prs_XPASS/XPASS_result/XPASS_AFR0815_rep1"
# pheno
pheno = read.table(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/',race,'/pheno_summary_out_GA/phenotypes_rho',rho,'_',GA,'.phen'),header=F,
                   colClasses = c(rep('character',2),rep('numeric',rep),rep('NULL',100-rep)),
                   skip = 110000)
pheno = as.data.frame(pheno[,c(1,rep+2)])
colnames(pheno) = c('id','y')

#####
# R^2 from XPASS
XPASS_file_path = paste0(file_path, "_PRS.txt")
PRS = fread(XPASS_file_path)
PRS = PRS[10001:20000,]
model <- lm(pheno$y~PRS$PRS_XPASS1)
r2 <- summary(model)$r.square
print(paste0("XPASS R^2 (PRS calculated by XPASS): ", r2))

# R^2 from plink
plink_file_path = paste0(file_path, "_PRS_plink.sscore")
PRS = fread(plink_file_path)
PRS = PRS[10001:20000,]
model <- lm(pheno$y~PRS$SCORE3_SUM)
r2 <- summary(model)$r.square
print(paste0("XPASS R^2 (PRS calculated by plink): ", r2))


