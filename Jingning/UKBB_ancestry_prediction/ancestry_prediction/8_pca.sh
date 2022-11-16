#!/usr/bin/env bash
#$ -N pca
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

module load python/3.9.10

python3 8_pca.py

