source('liklweitz_ghk_1.R')
surface_test = function(cost) {
  params = c(0.93, 0.66, cost, 0.26, -3)
  data = read.csv('genWeitzDataS45.csv', header = FALSE)
  D = 100
  seed = 45
  return(liklWeitz_ghk_1(params, data, D, seed))
}

log_c_grid = seq(0, 1, by = 0.01)
llk_values_R = sapply(log_c_grid, surface_test)

plot(log_c_grid, llk_values_R, type = "l", xlab = "log(c)", ylab = "Negative Log-Likelihood",
     main = "Likelihood Surface for Search Cost (R)")
abline(v = -3.13, col = "red", lty = 2)  # R’s estimate
abline(v = -3.01, col = "blue", lty = 2)  # MATLAB’s estimate