.libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages/",.libPaths())) #where winnerscurse installed

library(data.table)
setwd("/data/BB_Bioinformatics/Kevin/Confluence_GWAS/code")

#peaks are the independent loci we detected
peakfile="/vf/users/DCEG_Confluence/Kevin/finemap/detect_novel/Apr4.novel_loci_classified.tsv"
peaks=fread(peakfile)
outdir="/data/BB_Bioinformatics/HZ/confluence_result/"
#ancestry specific meta-analyses results:
eurmetalfile="../result/internal_ancestry_Feb10/Meta_Mar23_EUR1.tbl"
easmetalfile="../result/internal_ancestry_Feb10/Meta_Mar23_EAS1.tbl"
afrmetalfile="../result/internal_ancestry_Feb10/Meta_Mar23_AFR1.tbl"
amrmetalfile="../result/internal_ancestry_Feb10/Meta_Mar23_AMR1.tbl"
#to get the corrected beta
winners_curse=function(metalfile=eurmetalfile,pop="EUR")
{
  # #use winners_curse from Haoyu
  # #https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html
  library(winnerscurse)
  #winner's curse correction 1
  
  summary_data=as.data.frame(fread(metalfile))
  summary_data$rsid=summary_data$MarkerName
  summary_data$beta=summary_data$Effect
  summary_data$se=summary_data$StdErr
  summary_data=summary_data[summary_data$beta!=0 & summary_data$se!=0,]
  #it should only include "rsid"  "beta"  "se"
  summary_data1=summary_data[summary_data$MarkerName %in% peaks$MarkerName,]
  #only for all the independent variants
  out_CL <- conditional_likelihood(summary_data = summary_data1[,c("rsid","beta","se")],alpha = 1)
  save(out_CL,file=paste0(outdir,pop,"_out_CL.RData"))
  #out_CL = out_CL %>% select(rsid, beta.cl3)
  head(out_CL)
  
  #winner's curse correction 2
  #apply to full summary statistics
  out_FIQT <- FDR_IQT(summary_data = summary_data)
  save(out_FIQT,file=paste0(outdir,pop,"_out_FIQT.RData"))
  head(out_FIQT)
}
winners_curse(metalfile=eurmetalfile,pop="EUR")
winners_curse(metalfile=easmetalfile,pop="EAS")
winners_curse(metalfile=afrmetalfile,pop="AFR")
winners_curse(metalfile=amrmetalfile,pop="AMR")

compute_h2_winners=function(metalfile=eurmetalfile,opt="CL",pop="EUR")
{
  metalres=as.data.frame(fread(metalfile))
  peaks$SNP=peaks$MarkerName
  print(table(peaks$SNP %in% metalres$MarkerName))
  mypeaks=peaks[peaks$SNP %in% metalres$MarkerName,]
  idx=match(mypeaks$SNP,metalres$MarkerName)
  mymetalres=metalres[idx,]
  mymetalres$maf=mymetalres$Freq1
  mymetalres$maf[which(mymetalres$maf>0.5)]=1-mymetalres$maf[which(mymetalres$maf>0.5)]
  idx=mypeaks$novel
  knowndat=mymetalres[!idx,]
  noveldat=mymetalres[idx,]
  
  if (opt=="FIQT")
  {
    load(paste0(outdir,pop,"_out_FIQT.RData"))
    table(noveldat$MarkerName %in% out_FIQT$MarkerName)
    idx=match(noveldat$MarkerName,out_FIQT$MarkerName)
    if (sum(is.na(idx))>0)stop()
    noveldat$beta=out_FIQT$beta_FIQT[idx]
    table(knowndat$MarkerName %in% out_FIQT$MarkerName)
    idx=match(knowndat$MarkerName,out_FIQT$MarkerName)
    if (sum(is.na(idx))>0)stop()
    knowndat$beta=out_FIQT$beta_FIQT[idx]
  }
  
  if (opt=="CL")
  {
    load(paste0(outdir,pop,"_out_CL.RData"))
    table(noveldat$MarkerName %in% out_CL$rsid)
    idx=match(noveldat$MarkerName,out_CL$rsid)
    if (sum(is.na(idx))>0)stop()
    noveldat$beta=out_CL$beta.cl3[idx]
    table(knowndat$MarkerName %in% out_CL$rsid)
    idx=match(knowndat$MarkerName,out_CL$rsid)
    if (sum(is.na(idx))>0)stop()
    knowndat$beta=out_CL$beta.cl3[idx]
  }
  #2*f*(1-f)*(beta^2)
  novelh2=sum((noveldat$beta^2)*2*noveldat$maf*(1-noveldat$maf),na.rm=T)
  knownh2=sum((knowndat$beta^2)*2*knowndat$maf*(1-knowndat$maf),na.rm=T)
  h2=novelh2+knownh2
  return(data.frame(novelh2=novelh2,knownh2=knownh2,h2=h2))
}


get_res_winners=function(opt="FIQT")
{ 
  eur_h2=compute_h2_winners(metalfile = eurmetalfile,opt=opt,pop="EUR")
  amr_h2=compute_h2_winners(metalfile = amrmetalfile,opt=opt,pop="AMR")
  eas_h2=compute_h2_winners(metalfile = easmetalfile,opt=opt,pop="EAS")
  afr_h2=compute_h2_winners(metalfile = afrmetalfile,opt=opt,pop="AFR")
  
  tmp=data.frame(population=c("AFR","AMR","EAS","EUR"))
  tmp=cbind(tmp,rbind(afr_h2,amr_h2,eas_h2,eur_h2))
  #save output to file
  write.csv(tmp,file=paste0(outdir,"heritability_estimate_new_old_",opt,"_Apr10_useancestrymetalres.csv"),row.names = F,quote=F)
  
  return(tmp)
}

FIQT_res=get_res_winners(opt="FIQT")

result = read.csv(paste0(outdir,"heritability_estimate_new_old_",opt,"_Apr10_useancestrymetalres.csv"))

CL_res=get_res_winners(opt="CL")


compare_raw_corrected <- function(metalfile, pop) {
  metalres <- as.data.frame(fread(metalfile))
  peaks_sub <- peaks[peaks$MarkerName %in% metalres$MarkerName, ]
  idx <- match(peaks_sub$MarkerName, metalres$MarkerName)
  metal_sub <- metalres[idx, c("MarkerName", "Effect", "StdErr", "P-value", "Freq1")]
  
  load(paste0(outdir, pop, "_out_FIQT.RData"))
  
  idx2 <- match(metal_sub$MarkerName, out_FIQT$MarkerName)
  metal_sub$beta_FIQT <- out_FIQT$beta_FIQT[idx2]
  
  cat("==", pop, "==\n")
  cat("Matched peaks:", nrow(metal_sub), "\n")
  cat("Median abs raw beta:", median(abs(metal_sub$Effect), na.rm = TRUE), "\n")
  cat("Median abs corrected beta:", median(abs(metal_sub$beta_FIQT), na.rm = TRUE), "\n")
  cat("Mean abs raw beta:", mean(abs(metal_sub$Effect), na.rm = TRUE), "\n")
  cat("Mean abs corrected beta:", mean(abs(metal_sub$beta_FIQT), na.rm = TRUE), "\n")
  cat("Median SE:", median(metal_sub$StdErr, na.rm = TRUE), "\n")
  cat("Median -log10P:", median(-log10(metal_sub$"P-value"), na.rm = TRUE), "\n\n")
  
  return(metal_sub)
}

afr_sub <- compare_raw_corrected(afrmetalfile, "AFR")
eur_sub <- compare_raw_corrected(eurmetalfile, "EUR")
