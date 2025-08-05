simWeitz = function(N_cons, N_prod, param,  seed) {

  set.seed(seed)
  N_obs = N_cons * N_prod
  consumer = rep(1:N_cons, each = N_prod)
  N_prod_cons = rep(N_prod, N_cons)
  prod = rep(1:N_prod, N_cons)
  
  # Outside option (represents option of not buying) and brand fixed effects
  outside = (prod == 1)
  brand1 = (prod == 2)
  brand2 = (prod == 3)
  brand3 = (prod == 4)
  brand4 = (prod == 5)
  
  # Product characteristics: in this case only brand intercepts
  X = cbind(brand1, brand2, brand3, brand4)
  # First and last entry indices for each consumer
  index_first = match(unique(consumer), consumer)
  index_last = cumsum(N_prod_cons)
  # Parameters
  c = exp(param[length(param)]) * rep(1, N_obs)  # Search cost
  xb = rowSums(X * param[1:(length(param)-1)])  # Utility from observables
  # Draws affecting utility
  epsilon = rnorm(N_obs)
  eta = rnorm(N_obs)
  # Expected utility and utility
  eut = (xb + eta) * (1 - outside)
  ut = eut + epsilon
  
  
  # Forming z's:
  # Load lookup table
  table = read.csv("tableZ.csv", header = FALSE)
  # Initialize m
  m = numeric(N_obs)
  # Lookup table method - we omit the others for now
  for (i in 1:N_obs) {
    lookupvalue = abs(table[, 2] - c[i])
    if (table[1, 2] >= c[i] && c[i] >= table[nrow(table), 2]) {
      index_m = which.min(lookupvalue)
      m[i] = table[index_m, 1]
    } else if (table[1, 2] < c[i]) {
      m[i] = -c[i]
    } else if (c[i] < table[nrow(table), 2]) {
      m[i] = 4.001
    }
  }
  
  # Compute z
  z <- m + eut
  
  
  # Plug in large value for outside option as it is always "searched" first
  z = 100000 * outside + z * (1 - outside)
  
  # Create dataset
  da = data.frame(consumer = consumer, prod = prod, outside = outside, X = X, eut = eut, ut = ut, z = z)
  
  # Identify column indices
  whatz = ncol(da)
  whatu = whatz - 1
  whateu = whatu - 1
  
  # Sort data by z for each consumer
  values = numeric(N_obs)
  order = numeric(N_obs)
  for (i in 1:N_cons) {
    idx = index_first[i]:index_last[i]
    sorted = sort(da[idx, whatz], decreasing = TRUE, index.return = TRUE)
    values[idx] = sorted$x
    order[idx] = sorted$ix
  }
  
  # Compute reordering indices
  long_first = rep(index_first, N_prod_cons)
  order2 = order + long_first - 1
  
  # Reorder dataset
  data = da[order2, ]
  
  
  # Initialize search decision: 1. outside option always searched
  searched = outside
  
  # Search decision: 2. Search if z > all ut searched so far
  for (i in 1:N_cons) {
    for (j in (index_first[i] + 1):index_last[i]) {
      relevant_ut_sofar = data[index_first[i]:(j - 1), whatu] * searched[index_first[i]:(j - 1)]
      relevant_ut_sofar = relevant_ut_sofar[relevant_ut_sofar != 0]
      max_ut_sofar = if (length(relevant_ut_sofar) > 0) max(relevant_ut_sofar) else -Inf
      
      if (data[j, whatz] > max_ut_sofar) {
        searched[j] = 1
      }
    }
  }
  
  # Transaction: among those searched, pick max ut
  tran = rep(0, N_obs)
  searched_ut = data[, whatu] * searched
  
  for (i in 1:N_cons) {
    A = searched_ut[index_first[i]:index_last[i]]
    A[A == 0] = -100000
    indexch = which.max(A)
    tran[index_first[i] + indexch - 1] = 1
  }
  
  # Export data
  length = rep(N_prod, N_obs)  # Number of products per consumer
  searched_mat = matrix(searched, nrow = N_prod, ncol = N_cons)
  has_searched = searched_mat[2, ]  # Did consumer search at least once
  has_searched = matrix(rep(has_searched, each = N_prod), ncol = 1)
  last = rep(c(rep(0, N_prod - 1), 1), N_cons)  # Last product per consumer
  
  # Construct simulation matrix
  simulation = cbind(data[, 1:(whateu - 1)], last, has_searched, length, searched, tran)
  
  # Export to CSV
  filename = sprintf("genWeitzDataS%d.csv", seed)
  write.csv(simulation, filename, row.names = FALSE)
  
  
  return(simulation)
}