

## find the code of interested phenotype from ukbb html


## ancestry
a <- read.table("/dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab", nrows=1)
b <- as.character(a[1,])
which(b == "f.21000.0.0")
which(b == "f.22006.0.0")
cat /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab | awk -F '\t' '{ print $1, $9762, $9763, $9764, $9945}' > '/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/ethnic_background_from_ukbb.txt'


## HDL
a <- read.table("/dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab", nrows=1)
b <- as.character(a[1,])
which(b == "f.30760.0.0")
which(b == "f.30760.1.0")
mkdir /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/HDL
cat /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab | awk -F '\t' '{ print $1, $14456, $14457}' > '/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/HDL/HDL_from_ukbb.txt'

## LDL
a <- read.table("/dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab", nrows=1)
b <- as.character(a[1,])
which(b == "f.30780.0.0")
which(b == "f.30780.1.0")
mkdir /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/LDL
cat /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab | awk -F '\t' '{ print $1, $14484, $14485}' > '/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/LDL/LDL_from_ukbb.txt'



/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/EUR/chr1 \
--extract /users/jzhang2/tmp_tmp/snp \
--rm-dup exclude-all \
--threads 1 \
--make-bed \
--out /users/jzhang2/tmp_tmp/bfile

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--threads 1 \
--bfile /users/jzhang2/tmp_tmp/bfile \
--r bin4 \
--out /users/jzhang2/tmp_tmp/ld


/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--threads 1 \
--bfile /users/jzhang2/tmp_tmp/bfile \
--freq \
--out /users/jzhang2/tmp_tmp/freq
