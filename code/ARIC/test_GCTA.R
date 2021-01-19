setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/gcta_1.93.2beta")
system("gcta64 --bfile test --autosome --maf 0.01 --make-grm --out test --thread-num 10")