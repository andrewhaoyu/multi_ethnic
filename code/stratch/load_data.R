load("./data/missense.Rdata")

colnames(result)[90] = "STAAR"
rank_idx = order(as.numeric(result$STAAR))
result_new = result[rank_idx[1:3],]
as.numeric(result[,90])[971]
class(result)
result[971,90]
result$STAAR[971]
result[971,90]
