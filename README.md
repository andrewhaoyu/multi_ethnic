## Multi Ancestry PRS development
Simulation and Data Analyses for Multi-Ancestry PRS Development in Zhang et al. "Novel Methods for Multi-Ancestry Polygenic Prediction and their Evaluations in 5.1 Million Individuals of Diverse Ancestry." Biorxiv, 2022.

- Goal 

Polygenic risk scores (PRS) are becoming increasingly predictive of complex traits. However, their suboptimal performance in non-European ancestry populations raises questions about their clinical applications and the impact on health inequities. The goal of this study is to:
1. Develop a powerful and scalable method, CT-SLEB, based on ancestry-specific GWAS summary statistics from multi-ancestry training samples. The method integrates multiple techniques, including clumping and thresholding, empirical Bayes, and super learning.
2. Evaluate cutting-edge PRS approaches, including LDPred2, PRS-CSx, S-PolyPred+, etc.

All analyses are conducted by splitting the data into training, tuning, and validation datasets. The GWAS summary statistics from the training dataset are used to train the model. The tuning dataset is used to find the optimal tuning parameters. The final results are reported based on the independent validation dataset.


- Data 
1. The simulated data includes 600,000 subjects from five ancestries (AFR, AMR, EAS, EUR, SAS) generated using HapGen2 with reference data from the 1000 Genomes project. Each ancestry contains 120,000 subjects. The simulated genotype, phenotype, and GWAS summary statistics are available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP&version=4.1)
2. 23andMe data analyses. GWAS summary statistics for seven traits (any cardiovascular disease, depression, migraine diagnosis, singing back musical note, morning person, height, and heart metabolic disease burden) were obtained through collaboration with the research team at 23andMe.
3. The Global Lipids Genetics Consortium (GLGC) data analyses. GWAS summary statistics for four blood lipid traits are available at this [link](http://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/) (Graham, et al. Nature, 2021). 
4. The All of Us (AOU) data analyses. Two traits (height and BMI) were evaluated using data from All of Us.


- Code 
1. Code for simulation analyses can be found in three directories: `code/LD_simulation_GA`, `code/LD_simulation_large`, and `code/LD_simulation`. `code/LD_simulation_GA` contains the most up-to-date simulation analyses results, while the other two directories are more exploratory in nature as the project progresses.
2. Code for 23andMe data analyses can be found in `code/23andMe`.
3. Code for GLGC data analyses can be found in `code/GLGC_analysis`.
4. Code for AOU data analyses can be found in `code/AOU`.



