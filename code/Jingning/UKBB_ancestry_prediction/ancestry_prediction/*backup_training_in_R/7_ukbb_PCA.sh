#!/usr/bin/env bash
#$ -N ukbbpca
#$ -cwd
#$ -l mem_free=70G,h_vmem=70G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

R CMD BATCH --no-save --no-restore 7_ukbb_PCA.R
