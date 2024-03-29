load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/temp.rdata")
library(ggplot2)

source("/Users/zhangh24/GoogleDrive/multi_ethnic/code/LD_simulation_large/theme_Publication.R")
idx_sample = seq(1,nrow(prs_com),10)
data = prs_com[idx_sample,]

plot(data$effect_hm3,data$effect_mega)
abline(a = 0, b = 1, col = "red")
abline(h = 0, col = "red", lty = 2)
abline(v = 0, col = "red", lty = 2)


range(prs_com$effect_mega)
range(prs_com$effect_hm3)
idx_rank = order(abs(prs_com$effect_hm3),decreasing = T)
prs_com_order = prs_com[idx_rank,]
head(prs_com_order)
model =lm(effect_mega~effect_hm3,data = prs_com)
summary(model)
