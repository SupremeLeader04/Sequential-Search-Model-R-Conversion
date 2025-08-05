d = list.files(pattern = "^rezSimWeitz_ghk_D100S.*\\.csv$")
cumdata = c(rep(0, 10))
for (file in d) {
  data = head(read.csv(file, header=TRUE),-2)
  cumdata = cumdata + data[[1]]
}
avg = cumdata / 50
write.table(avg, "overall_results.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
