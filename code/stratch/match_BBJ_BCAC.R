#process the EAS summary statistics


library(data.table)
library(dplyr)
BCAC_Asian = as.data.frame(fread("BCAC_Asian/all.meta.icogs.oncoarray.asian.overall.Apr2016.1.txt"))
head(BCAC_Asian)

#create a unique identifier for alleles
#make sure the smaller character in the front
BCAC_Asian_update = BCAC_Asian %>% 
  mutate(unique_allele1 = toupper(ifelse(Allele1 < Allele2, Allele1,Allele2)),
         unique_allele2 = toupper(ifelse(Allele1 < Allele2, Allele2,Allele1)),
         chr_pos = sub("^([^_]+_[^_]+)_.*", "\\1", MarkerName),
         unique_SNP_id = paste0(chr_pos,"_",unique_allele1,"_",unique_allele2)) %>% 
  mutate(effect_allele = toupper(Allele1),
         non_effect_allele = toupper(Allele2),
         Freq_effect = Freq1,
         BETA = Effect,
         SE = StdErr,
         P = `P-value`,
         N_eff = 1/(2*Freq_effect*(1-Freq_effect)*SE^2)) %>% 
  select(unique_SNP_id, effect_allele, non_effect_allele, Freq_effect, BETA, SE, P, N_eff)

BBJ_Asian = as.data.frame(fread("East_Asian_GWAS_BBJ/hum0197.v3.BBJ.BC.v1/GWASsummary_BrC_Japanese_SakaueKanai2020.auto.txt"))
head(BBJ_Asian)

#clean BBJ Asian data to create an unique SNP ID to match with BCAC Asian data
BBJ_Asian_update = BBJ_Asian %>% 
  mutate(unique_allele1 = ifelse(Allele1 < Allele2, Allele1,Allele2),
         unique_allele2 = ifelse(Allele1 < Allele2, Allele2,Allele1),
         chr_pos = sub("^([^:]+):([^:]+).*", "\\1_\\2", v),
         unique_SNP_id = paste0(chr_pos,"_",unique_allele1,"_",unique_allele2)) %>% 
  mutate(effect_allele = toupper(Allele2),
         non_effect_allele = toupper(Allele1),
         #for effect sample size calculation, using the Allele frequency in controls
         Freq_effect = AF.Controls,
         P = p.value,
         N_eff = as.integer(1/(2*Freq_effect*(1-Freq_effect)*SE^2))) %>% 
  select(unique_SNP_id, effect_allele, non_effect_allele, Freq_effect, BETA, SE, P, N_eff, SNPID)
head(BBJ_Asian_update)


#rename the variable to prepare the merge step
BBJ_Asian_update_merge = BBJ_Asian_update  %>% 
  rename(effect_allele_BBJ = effect_allele, 
         non_effect_allele_BBJ = non_effect_allele, 
         Freq_effect_BBJ = Freq_effect, 
         BETA_BBJ = BETA, 
         SE_BBJ = SE,
         P_BBJ = P, 
         N_eff_BBJ = N_eff)  %>% 
  select(unique_SNP_id, effect_allele_BBJ, non_effect_allele_BBJ, Freq_effect_BBJ, BETA_BBJ, SE_BBJ, P_BBJ, N_eff_BBJ, SNPID)
BCAC_Asian_update_merge = BCAC_Asian_update  %>% 
  rename(effect_allele_BCAC = effect_allele, 
         non_effect_allele_BCAC = non_effect_allele, 
         Freq_effect_BCAC = Freq_effect, 
         BETA_BCAC = BETA, 
         SE_BCAC = SE,
         P_BCAC = P, 
         N_eff_BCAC = N_eff)  %>% 
  select(unique_SNP_id, effect_allele_BCAC, non_effect_allele_BCAC, Freq_effect_BCAC, BETA_BCAC, SE_BCAC, P_BCAC, N_eff_BCAC)
head(BBJ_Asian_update_merge)
head(BCAC_Asian_update_merge)

dim(BBJ_Asian_update_merge)
dim(BCAC_Asian_update_merge)
combine_data = full_join(BBJ_Asian_update_merge,BCAC_Asian_update_merge, by = "unique_SNP_id")
head(combine_data)

combine_data_update <- combine_data %>%
  mutate(
    # Check if both alleles exist and are different, then decide whether to flip the sign
    flip_sign = !is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) & (effect_allele_BBJ != effect_allele_BCAC),
    BETA_BCAC = ifelse(flip_sign, -BETA_BCAC, BETA_BCAC),
    # Proceed with the meta-analysis
    w_BBJ = ifelse(!is.na(BETA_BBJ), 1/(SE_BBJ^2), NA),
    w_BCAC = ifelse(!is.na(BETA_BCAC), 1/(SE_BCAC^2), NA),
    BETA_meta = ifelse(!is.na(w_BBJ) & !is.na(w_BCAC), 
                       (w_BBJ * BETA_BBJ + w_BCAC * BETA_BCAC) / (w_BBJ + w_BCAC), 
                       ifelse(is.na(w_BBJ), BETA_BCAC, BETA_BBJ)),
    SE_meta = ifelse(!is.na(w_BBJ) & !is.na(w_BCAC),
                     sqrt(1 / (w_BBJ + w_BCAC)),
                     ifelse(is.na(w_BBJ), SE_BCAC, SE_BBJ)),
    P_meta = 2 * pnorm(-abs(BETA_meta / SE_meta), lower.tail = TRUE)
  ) %>%
  mutate(
    # Set effect and non-effect allele for the meta-analysis
    effect_allele_meta = ifelse(!is.na(effect_allele_BBJ), effect_allele_BBJ, effect_allele_BCAC),
    non_effect_allele_meta = ifelse(!is.na(non_effect_allele_BBJ), non_effect_allele_BBJ, non_effect_allele_BCAC),
    # Calculate N_eff_meta
    N_eff_meta = ifelse(!is.na(N_eff_BBJ) & !is.na(N_eff_BCAC), 
                        N_eff_BBJ + N_eff_BCAC, 
                        ifelse(is.na(N_eff_BBJ), N_eff_BCAC, N_eff_BBJ))
  )
head(combine_data_update)
tail(combine_data_update)
write.table(combine_data_update, file = "Clean_summary_data/EAS_BCAC_BBJ_meta_sumdata.txt", row.names = F, col.names = T, quote = F)