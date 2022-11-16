#!/usr/bin/env bash
#$ -N pca5
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

module load python/3.9.10

python3 8_pca_5.py

