library(qrcode)
url <- "https://www.nature.com/articles/s41588-023-01501-z"
qr <- qr_code(url)
plot(qr)
png(filename = "ctsleb_qr.png")
plot(qr)
dev.off()
