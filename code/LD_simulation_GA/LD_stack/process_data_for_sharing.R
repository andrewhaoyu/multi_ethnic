## goal process genotype data for sharing ##
#load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.rdata")
# dup.id.list = list()
# library(data.table)
# for(i in 1:5){
#   data = fread(paste0(cur.dir,eth[i],"/duplicated.id"),header=F)
#   dup.id.list[[i]] = data
#   
# }
# 
# dup.id = unique(rbindlist(dup.id.list))
# colnames(dup.id) = "SNP"
#write.table(dup.id,file = paste0(cur.dir,"dup.id"),col.names = T,row.names = F,quote=F)
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
# args = commandArgs(trailingOnly = T)
# i = as.numeric(args[[1]])
# j = as.numeric(args[[2]])
# eth = c("EUR","AFR","AMR","EAS","SAS")
# 
# 
# library(data.table)
# library(dplyr)
# sid<-Sys.getenv('SLURM_JOB_ID')
# dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
# temp.dir = paste0('/lscratch/',sid,'/test/')
# cur.dir = "q"
# out.dir <-  paste0(cur.dir,eth[i],"/",eth[i],"_mega/")
# 
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.bed ",temp.dir,eth[i],"chr",j,".mega.bed"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.bim ",temp.dir,eth[i],"chr",j,".mega.bim"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.fam ",temp.dir,eth[i],"chr",j,".mega.fam"))
# 
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,"train.mega.bed ",temp.dir,eth[i],"chrtrain",j,".mega.bed"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,"train.mega.bim ",temp.dir,eth[i],"chrtrain",j,".mega.bim"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,"train.mega.fam ",temp.dir,eth[i],"chrtrain",j,".mega.fam"))
# 
# all_my_files <- rep("c",1)
#   all_my_files[1] = paste0(temp.dir,eth[i],"chrtrain",j,".mega")
# write.table(all_my_files,file = paste0(temp.dir,"all_my_files.txt"),row.names = F,col.names = F,quote=F)
# 
# res = system(paste0("/data/zhangh24/software/plink2 ",
#                     "--bfile ",temp.dir,eth[i],"chr",j,".mega ",
#                     "--merge-list ",temp.dir,"all_my_files.txt ",
#                     "--out ",temp.dir,"chr",j," ",
#                     "--exclude ",cur.dir,"dup.id ",
#                     "--make-bed "))
# 
# system(paste0("cd ",temp.dir,"; zip -r chr",j,".bed.zip chr",j,".bed"))
# system(paste0("cd ",temp.dir,"; zip -r chr",j,".bim.zip chr",j,".bim"))
# system(paste0("cd ",temp.dir,"; zip -r chr",j,".fam.zip chr",j,".fam"))
# system(paste0("cp ",temp.dir,"chr",j,".bed.zip ",out.dir,"chr",j,".bed.zip"))
# system(paste0("cp ",temp.dir,"chr",j,".fam.zip ",out.dir,"chr",j,".fam.zip"))
# system(paste0("cp ",temp.dir,"chr",j,".bim.zip ",out.dir,"chr",j,".bim.zip"))


# args = commandArgs(trailingOnly = T)
# i = as.numeric(args[[1]])
#eth = c("EUR","AFR","AMR","EAS","SAS")
# library(data.table)
# library(dplyr)
# cur.dir = "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# out.dir <-  paste0(cur.dir,eth[i],"/",eth[i],"_mega/")
# 
# system(paste0("cd " ,cur.dir,eth[i],"; ",
# "zip -r ",eth[i],"_mega.zip ",eth[i],"_mega"))
# 
# for(j in 1:22){
#   system(paste0("cd ",cur.dir,eth[i],"; ",
#                 "rm chr",j,"train.mega.* ; ",
#                 "rm chr",j,".mega.* "))
# }








############################process summary level statistics#########################
# args = commandArgs(trailingOnly = T)
# i = as.numeric(args[[1]])
# l = as.numeric(args[[2]])
# m = as.numeric(args[[3]])
# i1 = as.numeric(args[[4]])
# library(data.table)
# library(dplyr)
# 
# out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# n.rep = 10
# # for(i in 1:5){
# #   for(l in 1:3){
# #     for(m in 1:4){
# 
# #summary GA_3:5 only has mega+hm3 snps
# #summary GA_1:2 has 1KG snps
# #for sharing purpose, we restrct to mega+hm3 snps
# summmary.data.ref = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",1,"_GA_",3)))
# snp.ref = summmary.data.ref[,"SNP",drop=F]
# idx <- which(duplicated(snp.ref))
# if(length(idx)!=0){
#   snp.ref = snp.ref[-idx,,drop=F]  
# }
# 
# 
#       effect_mat_list = list()
#       names_list = list()
#       for(i_rep in 1:n.rep){
#         summary.data <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
#         idx <- which(duplicated(summary.data$SNP))
#         if(length(idx)!=0){
#           summary.data = summary.data[-idx,]  
#         }
#         summary.data = summary.data %>%
#           mutate(SE = BETA/STAT)
#         summary.data.select = left_join(snp.ref,summary.data)
#         coeff_mat = summary.data.select %>%
#           select(BETA,SE,P)
#         effect_mat_list[[i_rep]] = coeff_mat
#         names_list[[i_rep]] = paste0(c("BETA_rep_","SE_rep_","P_rep_"),i_rep)
#       }
#       snp.infor = summary.data.select %>%
#         select(SNP,CHR,BP,A1)
#       effect_mat = bind_cols(effect_mat_list)
#       names = as.vector(sapply(names_list,paste0))
#       colnames(effect_mat) = names
#       summary.out = cbind(snp.infor,effect_mat)
#       write.table(summary.out,file = paste0(out.dir.sum,eth[i],"/summary_combine/summary_rho_",l,"_size_",m,"_GA_",i1),
#                   row.names = F,col.names = T,quote=F)
#       system(paste0("cd " ,out.dir.sum,eth[i],"/summary_combine/ ; ",
#       "zip summary_rho_",l,"_size_",m,"_GA_",i1,".zip ",
#       "summary_rho_",l,"_size_",m,"_GA_",i1," ;",
#       "rm summary_rho_",l,"_size_",m,"_GA_",i1))

# #compress summary statistics#
# out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# for(i in 1:5){
#   system(paste0("cd " ,out.dir.sum,eth[i]," ; ",
#                       "zip -r ",eth[i],"_summary_combine.zip ",
#                       "summary_combine ;"))
# }
# # for(i in 1:5){
#   system(paste0("cd " ,out.dir.sum,eth[i]," ; ",
#                 "rm -rf summary_combine;"))
# }
################################################################################################

################################compress phenotypes data##########################################
# args = commandArgs(trailingOnly = T)

# library(data.table)
# library(dplyr)
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
# n.rep = 10
# for(i in 1:5){
#   for(l in 1:3){
#     for(i1 in 1:5){
#       y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
#       y <- y[,1:(2+n.rep),drop=F]
#       colnames(y) = c("FID","IID",paste0("pheno_",c(1:10)))
#       write.table(y,file = paste0(out.dir.sum,eth[i],"/pheno_combine/pheno_rho_",l,"_GA_",i1),
#                   row.names = F,col.names = T,quote=F)
#       system(paste0("cd " ,out.dir.sum,eth[i],"/pheno_combine/ ; ",
#                     "zip pheno_rho_",l,"_GA_",i1,".zip ",
#                     "pheno_rho_",l,"_GA_",i1," ;",
#                     "rm pheno_rho_",l,"_GA_",i1))
#       
#     }
#   }
# }
# for(i in 1:5){
#   system(paste0("cd " ,out.dir.sum,eth[i]," ; ",
#                       "zip -r pheno_combine.zip ",
#                       "pheno_combine ;",
#                 "rm -rf pheno_combine"))
# }



#     }
#   }
# }
######################################################################
      
#####################snp information##################################
# setwd("/data/zhangh24/multi_ethnic/")
# load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
# snp.infor.match.update = snp.infor.match %>% 
#   rename(SNP=id,A0=a0,A1=a1) %>% 
#   rename(FREQ_A1_AFR = AFR,
#          FREQ_A1_AMR = AMR,
#          FREQ_A1_EAS = EAS,
#          FREQ_A1_EUR = EUR,
#          FREQ_A1_SAS = SAS,
#          FREQ_A1_ALL = ALL) %>% 
#   select(SNP,CHR,position,TYPE,FREQ_A1_AFR,FREQ_A1_AMR,
#          FREQ_A1_EAS,FREQ_A1_EUR,FREQ_A1_SAS,FREQ_A1_ALL,rs_id)
# write.table(snp.infor.match.update,file = paste0("./result/LD_simulation_new/snp_infor"),
#             row.names = F,col.names = T,quote=F)
# snp.infor = fread(paste0("./result/LD_simulation_new/snp_infor"))
# #MEGA and HAPMAP snp list
# cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
# #mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
# #colnames(mega.infor)[5] <- "rsid"
# colnames(mega.list) = "rs_id"
# snp.infor.mega = inner_join(mega.list,snp.infor,by="rs_id")
# write.table(snp.infor.mega,file = paste0("./result/LD_simulation_new/snp_infor_mega+hm3"),
#             row.names = F,col.names = T,quote=F)
# 
# hm3.list =   as.data.frame(fread(paste0(cur.dir,"hm3rsid.txt")))[,1,drop=F]
# colnames(hm3.list) = "rs_id"
# snp.infor.hm3 = inner_join(hm3.list,snp.infor,by="rs_id")
# write.table(snp.infor.hm3,file = paste0("./result/LD_simulation_new/snp_infor_hm3"),
#             row.names = F,col.names = T,quote=F)
# 
# 
# pheno = fread("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/EUR/pheno_combine/pheno_rho_1_GA_1")
# sum = fread(paste0(out.dir.sum,eth[i],"/summary_combine/summary_rho_",l,"_size_",m,"_GA_",i1))
##################################snp information#######################################################













######################rename the summary statistics to reflect the ancestry prefix_
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i1 = as.numeric(args[[4]])
library(data.table)
library(dplyr)

out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
n.rep = 10
# for(i in 1:5){
#   for(l in 1:3){
#     for(m in 1:4){
#       for(i1 in 1:5){
        for(i_rep in 1:n.rep){
          filename = paste0("summary_rho_",l,"_size_",m,"_GA_",i1)
          update_filename = paste0(eth[i],"_",filename)
          command = paste0("cd ",out.dir.sum,eth[i],"/summary_combine/; ",
                           "unzip ",filename,".zip; ",
                           "mv ",filename," ",update_filename,"; ",
                           "zip ",update_filename,".zip ",update_filename,"; ",
                           "rm ",filename,".zip; ",
                           "rm ",update_filename)
          system(command)
        }
#########################################################################################


system(paste0("mkdir /data/zhangh24/multi_ethnic/result/LD_simulation_GA/summary_data"))
for(i in 1:5){
  command = paste0("cd ",out.dir.sum,eth[i],"; ",
                   "cp -r summary_combine /data/zhangh24/multi_ethnic/result/LD_simulation_GA/summary_data/",eth[i],"_summary_combine")
  system(command)
}

#       }
#     }
#   }
# }
#   
###############################################################################



####################compress summary statistics#########################################
# out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# for(i in 1:5){
#   system(paste0("cd " ,out.dir.sum,eth[i]," ; ",
#                       "zip -r ",eth[i],"_summary_combine.zip ",
#                       "summary_combine ;"))
# }
#  for(i in 1:5){
#   system(paste0("cd " ,out.dir.sum,eth[i]," ; ",
#                 "rm -rf summary_combine;"))
# }
################################################################################################







######################rename the phenotype data to reflect the ancestry prefix_

library(data.table)
library(dplyr)

out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
n.rep = 10
 for(i in 1:5){
   for(l in 1:3){
#     for(m in 1:4){
       for(i1 in 1:5){

  filename = paste0("pheno_rho_",l,"_GA_",i1)
  update_filename = paste0(eth[i],"_",filename)
  command = paste0("cd ",out.dir.sum,eth[i],"/pheno_combine/; ",
                   "unzip ",filename,".zip; ",
                   "mv ",filename," ",update_filename,"; ",
                   "zip ",update_filename,".zip ",update_filename,"; ",
                   "rm ",filename,".zip; ",
                   "rm ",update_filename)
  system(command)

       }
     }
   }

system(paste0("mkdir /data/zhangh24/multi_ethnic/result/LD_simulation_GA/phenotypes"))
for(i in 1:5){
  command = paste0("cd ",out.dir.sum,eth[i],"; ",
                   "cp -r pheno_combine /data/zhangh24/multi_ethnic/result/LD_simulation_GA/phenotypes/",eth[i],"_pheno_combine")
  system(command)
}
command =  paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_GA/;",
                  "zip -r phenotypes.zip phenotypes")
system(command)
##########################################################################




