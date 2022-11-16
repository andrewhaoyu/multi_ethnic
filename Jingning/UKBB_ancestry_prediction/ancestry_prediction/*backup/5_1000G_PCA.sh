#!/usr/bin/env bash
#$ -N pca
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

R CMD BATCH --no-save --no-restore 5_1000G_PCA.R
