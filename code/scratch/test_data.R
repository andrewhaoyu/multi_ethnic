# Reading a gzipped text file directly into R
data <- read.table(gzfile("/data/BB_Bioinformatics/ProjectData/BioBank_Japan/gsr_file_gz_female_auto/hum0197.v10.Autoimmune.ALL.txt.gz"), header = TRUE, sep = "\t")
idx <- which(data$POS==32609965)
data[idx,]
