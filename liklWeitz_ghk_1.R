liklWeitz_ghk_1 = function(param, data, D, seed) {
  
  # Extract consumer IDs
  consumer = data[,1]
  N_cons = length(unique(consumer))
  
  # Extract N_prod (third-to-last column, equivalent to MATLAB's data(:,end-2))
  N_prod = data[, ncol(data)-2]
  Js = unique(N_prod)
  Num_J = length(Js)
  consumerData = matrix(0, nrow = N_cons, ncol = 2)
  consumer_num = 0
  
  # Construct likelihood for consumers with same number of searches
  source("liklWeitz_ghk_2.R")
  
  for (i in 1:Num_J) {
    nalt = Js[i]
    dat = data[N_prod == nalt, , drop = FALSE]
    N_obs = nrow(dat)
    uniCons = N_obs / nalt
    consid2 = matrix(dat[,1], nrow = nalt, ncol = uniCons, byrow = FALSE)
    
    consumerData[(consumer_num + 1):(consumer_num + uniCons), 1] = t(consid2[1, ])
    consumerData[(consumer_num + 1):(consumer_num + uniCons), 2] = liklWeitz_ghk_2(param, dat, D, nalt, seed)
    consumer_num = consumer_num + uniCons
  }
  
  # Sum over consumers
  llk = -sum(log(consumerData[,2]))
  
  # Check for errors or save output
  if (is.na(llk) || llk == Inf || llk == -Inf || !is.numeric(llk) || is.complex(llk)) {
    loglik = 1e+300
  } else {
    loglik = llk
    print(param)
    paramLL = c(param, loglik)
    # Save preliminary output
    write.csv(t(paramLL), sprintf("betaWeitz_ghk_D%dS%d.csv", D, seed), row.names = FALSE)
  }
  
  return(loglik)
}