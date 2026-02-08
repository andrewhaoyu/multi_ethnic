library(devtools)
install_github("amandaforde/winnerscurse")




library(winnerscurse)
## simulated GWAS summary statistics
sim_dataset <- sim_stats(nsnp=10^6,h2=0.4,prop_effect=0.01,nid=50000)
summary_stats <- sim_dataset$dissim_datasetsummary_stats <- sim_dataset$disc


#winner's curse correction 1
#apply to full summary statistics
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha = 5e-8) 
head(out_CL)
out_CL = out_CL %>% select(rsid, beta.cl3)
head(out_CL)

#winner's curse correction 2
#apply to full summary statistics
out_FIQT <- FDR_IQT(summary_data = summary_stats)
head(out_FIQT)

#generate allele frequency
summary_stats$f = runif(10^6,0,0.5)

#only apply to significant variants
#merge the data to significant variants
library(dplyr)
summary_stats_sig = left_join(out_CL,summary_stats,by = "rsid")
out_FIQT = out_FIQT %>% select(rsid,beta_FIQT)
summary_stats_sig = left_join(summary_stats_sig, out_FIQT)
beta = summary_stats_sig$beta
se = summary_stats_sig$se
f = summary_stats_sig$f

#option1: calculating variance correcting for sampling error
#E(X)^2 = E(X^2)-var(X)
#per variant variance = 2*f*(1-f)*(beta^2-se^2)

summary_stats_sig$herit = 2*f*(1-f)*(beta^2-se^2)

sum(summary_stats_sig$herit)

#option2: correct for winner's curse using conditional likelihood
beta_cl = out_CL$beta.cl3
summary_stats_sig$herit_cl = 2*f*(1-f)*(beta_cl^2)
sum(summary_stats_sig$herit_cl)

#option3: correct for winner's curse using FDR inverse quantile transformation
beta_FIQT = summary_stats_sig$beta_FIQT
summary_stats_sig$herit_FIQT = 2*f*(1-f)*(beta_FIQT^2)
sum(summary_stats_sig$herit_FIQT)
