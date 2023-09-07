###########################################################
# Individual analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(devtools)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
#install_github("zilinli1988/SCANG")
#install_github("xihaoli/STAARpipeline")
library(STAARpipeline)

#Get Input array
args=commandArgs(TRUE)

###########################################################
#           User Input
###########################################################
## Number of jobs for each chromosome
jobs_num <- 1
## aGDS directory
agds_dir <- "/data/COVID_ADHOC/lxwg/STAARpipeline/Steps/Association_Analysis/Step2/Input/agds_dir.Rdata"
## Null model
obj_nullmodel <- "/data/COVID_ADHOC/lxwg/STAARpipeline/Steps/Association_Analysis/Step2/Input/obj_nullmodel.Rdata"
## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
output_path <- args[4] #"/path_to_the_output_file/"
## output file name
output_file_name <- args[5] #"TOPMed_F5_LDL_Individual_Analysis"

## input array id from batch file (Harvard FAS RC cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)
load(agds.path)
start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc <= end_loc)
{
  results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                     genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=20,
                                                     QC_label=QC_label,variant_type=variant_type,
                                                     geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

