# module load python/3.9.10
# compute-054

# Two measurements at the initial visit (*_0_0 and *_0_1) are averaged
import pandas as pd
import os
import sys
# os.environ["PYSPARK_SUBMIT_ARGS"]="--total-executor-cores 10 --executor-memory 50g --driver-memory 50g pyspark-shell"
import pickle
import os

os.chdir("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/pca_5")
 
import hail as hl
hl.init()

## Path to input files
bed_path = "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall.bed"
bim_path = "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall.bim"
fam_path = "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall.fam"

## Import imputed genotype data
dataset = hl.import_plink(bed_path, bim_path, fam_path)

eigenvalues, scores, loadings = hl.hwe_normalized_pca(dataset.GT, k=5)

fileObj = open('eigenvalues.obj', 'wb')
pickle.dump(eigenvalues, fileObj)
fileObj.close()

fileObj = open('scores.obj', 'wb')
pickle.dump(scores, fileObj)
fileObj.close()

fileObj = open('loadings.obj', 'wb')
pickle.dump(loadings, fileObj)
fileObj.close()


#
#
# htt, rf_model = assign_population_pcs(
#     scores,
#     pc_cols = ht.scores,
#     known_col,
# )
#
#
#
#
#
# import pickle
# from gnomad.sample_qc.ancestry import assign_population_pcs
#
# # Load MatrixTable for projection and gnomAD loadings Hail Table
# loadings_ht = hl.read_table(path_to_gnomad_loadings)
#
# # Project new genotypes onto loadings
# ht = hl.experimental.pc_project(
#     mt_to_project.GT,
#     loadings_ht.loadings,
#     loadings_ht.pca_af,
# )
#
#
# # Assign global ancestry using the gnomAD RF model and PC project scores
# # Loading of the v2 RF model requires an older version of scikit-learn, this can be installed using pip install -U scikit-learn==0.21.3
# with hl.hadoop_open(path_to_gnomad_rf, "rb") as f:
#     fit = pickle.load(f)
#
# # Reduce the scores to only those used in the RF model, this was 6 for v2 and 16 for v3.1
# num_pcs = fit.n_features_
# ht = ht.annotate(scores=ht.scores[:num_pcs])
# htt, rf_model = assign_population_pcs(
#     ht,
#     pc_cols = ht.scores,
#     fit = fit,
# )
#
# # The returned Hail Table includes the imputed population labels and RF probabilities for each gnomAD global population
#
