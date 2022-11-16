#!/usr/bin/env bash
#$ -N allchr_bim
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

R CMD BATCH --no-save --no-restore bim.R

