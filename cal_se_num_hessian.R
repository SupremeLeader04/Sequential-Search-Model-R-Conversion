cal_se_num_hessian = function(param, step, data, D, seed) {
  
  # Extract consumer IDs
  consumer = data[,1]
  N_cons = length(unique(consumer))
  
  # Extract N_prod (assuming it's the third-to-last column, adjust index as needed)
  N_prod = data[, ncol(data)-2]
  Js = unique(N_prod)
  Num_J = length(Js)
  
  # Initialize vectors
  consumerData = numeric(N_cons)
  consumerData_new = numeric(N_cons)
  consumer_num = 0
  
  # Construct likelihood for consumers with the same number of searches
  source("liklWeitz_ghk_2.R")
  for (i in 1:Num_J) {
    nalt = Js[i]
    dat = data[N_prod == nalt,]
    N_obs = nrow(dat)
    uniCons = N_obs / nalt
    
    consumerData[(consumer_num+1):(consumer_num+uniCons)] = 
      liklWeitz_ghk_2(param, dat, D, nalt, seed)
    consumer_num = consumer_num + uniCons
  }
  
  # Handle small values and take log
  consumerData_t = consumerData < 4.9407e-324
  consumerData[consumerData_t] = consumerData[consumerData_t] + 4.9407e-324
  consumerData = log(consumerData)
  
  # Initialize gradient matrix
  grad_mat = matrix(0, nrow = N_cons, ncol = length(param))
  
  # Take steps away from each parameter and construct likelihood
  for (p in 1:length(param)) {
    param_step = param
    param_step[p] = param[p] + step
    consumer_num = 0
    for (i in 1:Num_J) {
      nalt = Js[i]
      dat = data[N_prod == nalt,]
      N_obs = nrow(dat)
      uniCons = N_obs / nalt
      
      consumerData_new[(consumer_num+1):(consumer_num+uniCons)] = 
        liklWeitz_ghk_2(param_step, dat, D, nalt, seed)
      consumer_num = consumer_num + uniCons
    }
    
    # Handle small values and take log
    consumerData_new_t = consumerData_new < 4.9407e-324
    consumerData_new[consumerData_new_t] = consumerData_new[consumerData_new_t] + 4.9407e-324
    consumerData_new = log(consumerData_new)
    
    # Compute gradient for parameter p
    grad_mat[, p] = (consumerData_new - consumerData) / step
  }
  
  # Compute information matrix and standard errors
  I = t(grad_mat) %*% grad_mat / N_cons
  se = sqrt(diag(solve(I * N_cons)))
  
  return(se)
}