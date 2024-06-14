#CT-SLEB for GLGC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#l represent trait

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
library(data.table)
library(dplyr)
library(devtools)
# cd /home/$USER/R/
#mkdir 4.4
#cd 4.4
#mkdir library
#export R_LIBS_USER=/home/$USER/R/4.4/library
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
install_github("andrewhaoyu/CTSLEB")
library(CTSLEB,lib.loc = "/home/zhangh24/R/4.4/library/")
eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
eth_EUR = "EUR"
eth = eth_vec[i]
trait = trait_vec[l]
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth_EUR,"/all_chr.bed ",temp.dir,eth_EUR,"all_chr.bed"))
system(paste0("cp ",kg.dir,eth_EUR,"/all_chr.bim ",temp.dir,eth_EUR,"all_chr.bim"))
system(paste0("cp ",kg.dir,eth_EUR,"/all_chr.fam ",temp.dir,eth_EUR,"all_chr.fam"))

system(paste0("cp ",kg.dir,eth,"/all_chr.bed ",temp.dir,eth,"all_chr.bed"))
system(paste0("cp ",kg.dir,eth,"/all_chr.bim ",temp.dir,eth,"all_chr.bim"))
system(paste0("cp ",kg.dir,eth,"/all_chr.fam ",temp.dir,eth,"all_chr.fam"))
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/Ruoyu/",eth,"/",trait,"/")

#load EUR data
sum_eur = as.data.frame(fread(paste0(data.dir,"EUR/",trait,".txt"),header=T))


sum_eur = sum_eur %>% 
  select(rsID, CHR, POS_b37, BETA, SE, A1, P) %>% 
  rename(SNP = rsID, BP = POS_b37)
#load target population data
sum_tar = as.data.frame(fread(paste0(data.dir,eth,"/",trait,".txt"),header=T))
sum_tar = sum_tar %>% 
  select(rsID, CHR, POS_b37, BETA, SE, A1, P) %>% 
  rename(SNP = rsID, BP = POS_b37)
#align allels
sum_com <- AlignSum(sum_target = sum_tar,
                    sum_ref = sum_eur)
SplitSum <- function(x,
                     results_dir,
                     write_tables = TRUE){
  print("executing SplitSum() ...")
  sum_com <- x
  sum_com_select <- sum_com %>%
    mutate(split_ind =
             ifelse(
               (P<P_ref)|is.na(P_ref),T,F)
    )%>%
    select(SNP,P,P_ref,split_ind)
  
  sum_com_select_other_ref <- sum_com_select %>%
    filter(split_ind==F) %>%
    select(SNP,P_ref) %>%
    rename(P = P_ref)
  
  sum_com_select_tar_ref <- sum_com_select %>%
    filter(split_ind==T) %>%
    select(SNP,P)
  
  sum_list <- list(sum_com_select_other_ref,
                   sum_com_select_tar_ref)
  
 
  return(sum_list)
  
}

#split the SNPs into two groups
sum_com_split <- SplitSum(sum_com)
#sum_com_split is a list with two data frame
#the first data frame contains SNPs with p_eur < p_target, the p-value column is from p_eur.
sum_other_ref = sum_com_split[[1]]
#the second data frame contains target population-specific SNPs or p_target < p_eur. The p-value column is from p_target. 
sum_tar_ref = sum_com_split[[2]]




###################CLUMPING STEP started##############

# we use plink1.9 for the clumping purpose
# specify vector for clumping r square and base window size
#the clumping_window_ize = base_window_size/clumping_r_square so that lower clumping r2 can have larger clumping window size
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

write.table(sum_other_ref,paste0(temp.dir,"sum_other_ref"),col.names = T,row.names = F,quote=F)
write.table(sum_tar_ref,paste0(temp.dir,"sum_tar_ref"),col.names = T,row.names = F,quote=F)
# --clump-p1 determines the upper bound of p-value to be kept in the clumping. We set it as 1.
# --clump-r2 is the clumping r square
# --clump-kb is the clumping window size
snp_list = list()
temp = 1
soft.dir = "/data/zhangh24/software/"
snp_list = list()
temp = 1
for(r_ind in 1:length(r2_vec)){
  #create the window size given the clumping r2
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    pthr = 1.0
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    #for the first group, we perform clumping using EUR popultion as the reference   
    system(paste0(soft.dir,"plink2 ",
                  "--bfile ",temp.dir,eth_EUR,"all_chr ",
                  "--clump ",temp.dir,"sum_other_ref ",
                  "--clump-p1 ",pthr," ",
                  "--clump-r2 ",r2thr," ",
                  "--clump-kb ",kbpthr," ",
                  "--out ", temp.dir,"EUR_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    #for the second group, we perform clumping using AFR population as the reference
    system(paste0(soft.dir,"plink2 ",
                  "--bfile ",temp.dir,eth,"all_chr ",
                  "--clump ",temp.dir,"sum_tar_ref ",
                  "--clump-p1 ",pthr," ",
                  "--clump-r2 ",r2thr," ",
                  "--clump-kb ",kbpthr," ",
                  "--out ", temp.dir,eth,"_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    
    #combine the SNPs from the two clumping groups
    LD_EUR= fread(paste0(temp.dir,"EUR_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD_tar = fread(paste0(temp.dir,eth,"_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD  = rbind(LD_EUR,LD_tar)
    snp_list[[temp]] = LD
    names(snp_list[[temp]]) = paste0("clump_r2_",r2thr,"_ws_",kbpthr)
    temp = temp + 1
  }
}

save(snp_list, file = paste0(out.dir,"snp_list.rdata"))