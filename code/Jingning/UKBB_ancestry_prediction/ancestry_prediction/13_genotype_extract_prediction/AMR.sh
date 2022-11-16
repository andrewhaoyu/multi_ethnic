#!/usr/bin/env bash
#$ -N ukbb_geno_AMR
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=5000G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/
mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data

ethnic='AMR'

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/${ethnic}

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr${SGE_TASK_ID} \
--keep /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/${ethnic}.all_id \
--extract /dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt \
--rm-dup exclude-all \
--threads 1 \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/${ethnic}/chr${SGE_TASK_ID}
