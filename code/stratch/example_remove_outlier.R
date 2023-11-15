data <- mtcars # Example dataset
mean_mpg <- mean(data$mpg)
sd_mpg <- sd(data$mpg)
cutoff <- 3 * sd_mpg # 3 standard deviations

# Identify outliers
outliers <- data$mpg > mean_mpg + cutoff | data$mpg < mean_mpg - cutoff

# Remove outliers
data_clean <- data[!outliers, ]
