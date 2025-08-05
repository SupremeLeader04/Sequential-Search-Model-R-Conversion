d = list.files(pattern = "^rezSimWeitz_ghk_D100S.*\\.csv$")
num_files = length(d)
data_matrix = matrix(0, nrow = 12, ncol = num_files)

for (i in 1:num_files) {
  file = d[i]
  data = read.csv(file, header = TRUE)[1:12, 1]
  data_matrix[, i] = data
}

M = numeric(10)  
M2 = numeric(10) 

for (i in 1:10) {
  M[i] = mean(data_matrix[i, ])
  print(data_matrix[i,])
  M2[i] = sd(data_matrix[i, ])
}

write.table(M, "result_ghk_mean_Sum.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M2, "result_ghk_std_Sum.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)