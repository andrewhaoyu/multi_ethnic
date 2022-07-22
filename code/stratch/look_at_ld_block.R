library(ggplot2)
library(data.table)
library(dplyr)
data1 = fread("./data/ld_block/ld_block_EUR")
data2 = fread("./data/ld_block/ld_block_AFR")
ggplot(data = data) + 
  geom_segment(aes(x = start, xend = stop, y = start, yend = stop))


