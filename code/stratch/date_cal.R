mydates <- as.Date(c("2022-04-05", "2014-09-01"))
# number of days between 6/22/07 and 2/13/04
days <- mydates[1] - mydates[2]
days = as.numeric(days)

ConvertYMD <- function(days){
  year = as.integer(days/365)
  month = as.integer((days - 365*year)/30)
  day = as.integer(days - year*365-month*30)
  print(c(year,month,day))
}
ConvertYMD(days)
