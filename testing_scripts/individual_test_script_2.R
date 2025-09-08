library(dplyr)

# SAFE READ CSV TO HANDLE EMPTY FILES LIKE MATLAB
safe_read_csv <- function(filename, header = FALSE) {
  if (file.size(filename) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))  # empty matrix
  } else {
    return(read.csv(filename, header = header))
  }
}

for (s in 1:1) {
  # liklweitz1 parts
  seed = s
  param = c(0.92752, 0.62911, 0.43017,0.2359, -3.0017)
  D = 100
  data = safe_read_csv(sprintf("genWeitzDataS%d.csv", seed), header = FALSE)
  
  # Extract consumer IDs
  consumer = data[,1]
  N_cons = length(unique(consumer))
  
  # Extract N_prod (third-to-last column, equivalent to MATLAB's data(:,end-2))
  N_prod = data[, ncol(data)-2]
  Js = unique(N_prod)
  Num_J = length(Js)
  consumerData = matrix(0, nrow = N_cons, ncol = 2)
  consumer_num = 0
  
  nalt = Js
  dat = data[N_prod == nalt, , drop = FALSE]
  N_obs = nrow(dat)
  uniCons = N_obs / nalt
  consid2 = matrix(dat[,1], nrow = nalt, ncol = uniCons, byrow = FALSE)
  consumerData[(consumer_num + 1):(consumer_num + uniCons), 1] = t(consid2[1, ])
  
  # liklweitz2 parts
  consumer = dat[,1]
  N_obs = length(consumer)
  N_cons = length(unique(consumer))
  N_prod = nalt
  
  tran = dat[,ncol(dat)]
  searched = dat[,ncol(dat)-1]
  tran_ranking = colSums(matrix(tran, nrow=N_prod, ncol=N_cons) * rep(1:N_prod, N_cons))
  searched_amt = colSums(matrix(searched, nrow=N_prod, ncol=N_cons))
  last_searched = matrix(0, nrow=N_prod, ncol=N_cons)
  for(i in 1:N_cons) {
    last_searched[searched_amt[i], i] = 1
  }
  
  c = exp(param[length(param)]) * rep(1, N_obs)
  X = dat[,4:(3+length(param[1:(length(param)-1)]))]
  xb = rowSums(X * matrix(param[1:4], nrow=nrow(X), ncol=4, byrow=TRUE))
  
  table = safe_read_csv('tableZ.csv', header=FALSE)
  m = rep(0, N_obs)
  for(i in 1:N_obs) {
    lookupvalue = abs(table[,2] - c[i])
    if(table[1,2] >= c[i] && c[i] >= table[nrow(table),2]) {
      index_m = which.min(lookupvalue)
      m[i] = table[index_m,1]
    } else if(table[1,2] < c[i]) {
      m[i] = -c[i]
    } else if(c[i] < table[nrow(table),2]) {
      m[i] = 4.001
    }
  }
  
  etaDraw = as.matrix(safe_read_csv(sprintf("etaDrawS%d.csv", seed), header=FALSE))
  xb_con_prod = matrix(xb, nrow=N_prod, ncol=N_cons)
  mc = matrix(m, nrow=N_prod, ncol=N_cons)
  L_ALL = matrix(0, nrow=N_cons, ncol=D)
  
  # detecting individual cases
  
  # case3_consumer = 0
  # 
  # for(con in 1:N_cons) {
  #   H = searched_amt[con]                               # The last searched product
  #   mu_i = matrix(0, nrow=N_prod, ncol=D)               # Presearch shock vector
  #   eps_i = matrix(0, nrow=N_prod, ncol=D)              # Postsearch shock vector
  #   S_bar = setdiff(1:N_prod, 1:H)                      # The set of unsearched products
  #   tr = tran_ranking[con]                              # Position of purchased product
  #   # Step 1: Sampling mu_i0
  #   mu_i[1, ] = etaDraw[con, ]  # Sampling for the outside option
  # 
  #   if(tr == 1 && H > 1) {
  #     case4_consumer = con
  #     print(case4_consumer)
  #     print(H)
  #   }
  # }
  
  # probs = c(74, 97, 116, 218, 295, 393, 414, 442, 496, 546, 708, 723, 729, 750, 841, 848, 883, 918, 939)
  
  case = numeric(1000)
  
  for (con in 1:1000) {
    
    H = searched_amt[con]                               # The last searched product
    mu_i = matrix(0, nrow=N_prod, ncol=D)               # Presearch shock vector
    eps_i = matrix(0, nrow=N_prod, ncol=D)              # Postsearch shock vector
    S_bar = setdiff(1:N_prod, 1:H)                      # The set of unsearched products
    tr = tran_ranking[con]                              # Position of purchased product
    
    
    # Step 1: Sampling mu_i0
    mu_i[1, ] = etaDraw[con, ]  # Sampling for the outside option
    
    
    # Case 1 (main paper): consumer searched products and only purchased the outside option
    if(tr == 1 && H > 1) {
      
      case[con] = 1
      
      # Step 2: Sampling mu_iH
      b_iH = mu_i[1, ] - xb_con_prod[H, con] - mc[H, con]
      randvals_30_1 = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_1.csv", seed, con), header = FALSE))
      mu_i[H, ] = qnorm(randvals_30_1 * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      if (H >= 3) {
        for(h in (H - 1):2) {
          b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
          randvals_30_2 = as.numeric(safe_read_csv(sprintf('ranvals_S%d_Con%d_2%d.csv', seed, con, h), header = FALSE))
          mu_i[h, ] = qnorm(randvals_30_2 * (1 - pnorm(b_ih)) + pnorm(b_ih))
        }
      }
      
      
      # Stopping rule
      l_stop = rep(1, D)
      for(l in S_bar) {
        l_stop = l_stop * pnorm(mu_i[1, ] - xb_con_prod[l, con] - mc[l, con])
      }
      
      # Continuation rule
      l_cont = 1 - pnorm(b_iH)
      
      # Selection rule
      l_selection = rep(1, D)
      if (H >= 3) {
        for(h in (H-1):2) {
          b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
          l_selection = l_selection * (1 - pnorm(b_ih))
        }
      }
      
      # Choice rule
      l_choice = rep(1, D)
      for(h in 2:H) {
        l_choice = l_choice * pnorm(mu_i[1, ] - mu_i[h, ] - xb_con_prod[h, con])
      }
    
  
    
    
    
    
  
    # Case 2 (Appendix A): consumer searched products and purchased the product searched last
    } else if(tr == H && H > 1) {
      
      case[con] = 2
      
      # Step 1-2: Sampling mu_il
      if (length(S_bar) > 0) {
        ranvals_3_1 = as.matrix(safe_read_csv(sprintf("ranvals_S%d_Con%d_1.csv", seed, con), header=FALSE))
        mu_i[S_bar, ] = ranvals_3_1
      }
      
      # Step 2: Sampling mu_iH
      upper_bound = rbind(mu_i[1, ] - xb_con_prod[H, con] - mc[H, con],
                          -xb_con_prod[H, con] + xb_con_prod[S_bar, con] + mu_i[S_bar, ])
      b_iH = apply(upper_bound, 2, max)
      ranvals_3_2 = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_2.csv", seed, con), header=FALSE))
      mu_i[H, ] = qnorm(ranvals_3_2 * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      if (H >= 3) {
        for(h in (H-1):2) {
          b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
          ranvals_3_3h = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_3%d.csv", seed, con, h), header=FALSE))
          mu_i[h, ] = qnorm(ranvals_3_3h * (1 - pnorm(b_ih)) + pnorm(b_ih))
        }
      }
      
      # Step 4: Sampling epsilon_i(1:H - 1)
      if (H >= 3) {
        for(h in 2:(H-1)) {
          b_eps_ih = xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[h, con] - mu_i[h, ]
          ranvals_3_4h = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_4%d.csv", seed, con, h), header=FALSE))
          eps_i[h, ] = qnorm(ranvals_3_4h * pnorm(b_eps_ih))
        }
      }
      
      
      # Stopping rule
      if (H > 2) {
        middle_term <- mu_i[2:(H-1), ] + xb_con_prod[2:(H-1), con] + eps_i[2:(H-1), ]
        l_stop <- 1 - pnorm(apply(rbind(mu_i[1, ],
                                        (xb_con_prod[S_bar, con] + mu_i[S_bar, ] + mc[S_bar, con]),
                                        middle_term),
                                  2, max) - xb_con_prod[H, con] - mu_i[H, ])
      } else {
        l_stop <- 1 - pnorm(apply(rbind(mu_i[1, ],
                                        (xb_con_prod[S_bar, con] + mu_i[S_bar, ] + mc[S_bar, con])),
                                  2, max) - xb_con_prod[H, con] - mu_i[H, ])
      }
      
      # Continuation rule
      l_cont = 1 - pnorm(b_iH)
      
      # Selection rule
      l_selection = rep(1, D)
      if (H >= 3) {
        for(h in seq(H-1, 2, by = -1)) {
          b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
          l_selection = l_selection * (1 - pnorm(b_ih))
        }
      }
      
      # Choice rule
      l_choice = rep(1, D)
      if (H >= 3) {
        for(h in 2:(H-1)) {
          l_choice = l_choice * pnorm(xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[h, con] - mu_i[h, ])
        }
      }
  
  
      
      
      
      
      # Case 3 (Appendix A): consumer searched products and purchased a product neither not searched last nor outside option
    } else if((tr > 1) && (tr < H)) {
      
      case[con] = 3
      
      
      # Step 1-2: Sampling mu_il
      ranvals_2_1 = as.matrix(safe_read_csv(sprintf("ranvals_S%d_Con%d_1.csv", seed, con), header=FALSE))
      mu_i[S_bar, ] = ranvals_2_1
      
      # Step 2: Sampling mu_iH
      upper_bound = rbind(mu_i[1, ] - xb_con_prod[H, con] - mc[H, con],
                          -xb_con_prod[H, con] + xb_con_prod[S_bar, con] + mu_i[S_bar, ])
      b_iH = apply(upper_bound, 2, max)
      ranvals_2_2 = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_2.csv", seed, con), header=FALSE))
      mu_i[H, ] = qnorm(ranvals_2_2 * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      for(h in (H-1):2) {
        b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
        ranvals_2_3h = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_3%d.csv", seed, con, h), header=FALSE))
        mu_i[h, ] = qnorm(ranvals_2_3h * (1 - pnorm(b_ih)) + pnorm(b_ih))
      }
      
      
      # Step 4: Sampling epsilon for the purchased option
      b_eps_upper = xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[tr, con] - mu_i[tr, ]
      max_values = rbind(xb_con_prod[S_bar, con] + mu_i[S_bar, ] + mc[S_bar, con],
                         mu_i[1, ])
      b_eps_lower = apply(max_values, 2, max) - xb_con_prod[tr, con] - mu_i[tr, ]
      ranvals_2_4 = as.numeric(safe_read_csv(sprintf("ranvals_S%d_Con%d_4.csv", seed, con), header=FALSE))
      eps_i[tr, ] = qnorm(ranvals_2_4 * (pnorm(b_eps_upper) - pnorm(b_eps_lower)) + pnorm(b_eps_lower))
      
      
      # Stopping rule
      l_stop = pnorm(b_eps_upper) - pnorm(b_eps_lower)
      
      # Continuation rule
      l_cont = 1 - pnorm(b_iH)
      
      
      # Selection rule
      l_selection = rep(1, D)
      for(h in (H-1):2) {
        b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
        l_selection = l_selection * (1 - pnorm(b_ih))
      }
      
      
      # Choice rule
      l_choice = rep(1, D)
      for(h in setdiff(2:H, tr)) {
        l_choice = l_choice * pnorm(xb_con_prod[tr, con] + mu_i[tr, ] + eps_i[tr, ] - xb_con_prod[h, con] - mu_i[h, ])
      }
      
      
      
      
      
      
      # Case 4 (Appendix A): consumer only searched outside option and purchased the outside option
    } else if(tr == 1 && H == 1) {
      
      case[con] = 4
      
      # Stopping rule
      l_stop = rep(1, D)
      for(l in S_bar) {
        l_stop = l_stop * pnorm(mu_i[1, ] - xb_con_prod[l, con] - mc[l, con])
      }
      
      # Continuation rule
      l_cont = rep(1, D)
      
      # Selection rule
      l_selection = rep(1, D)
      
      # Choice rule
      l_choice = rep(1, D)
    }
  
  
    L_ALL[con, ] = l_stop * l_cont * l_selection * l_choice
  }
  
  
  llk = rowMeans(L_ALL)
  
  # ### COMPARING ###
  # 
  # r_llk = log(llk)
  # mat_llk = read.csv(sprintf('mat_likl_S%d.csv', seed), header = F)[[1]]
  # 
  # analysis = data.frame(
  #   'con' = 1:N_cons,
  #   'case' = case,
  #   'r_llk' = r_llk,
  #   'mat_llk' = mat_llk,
  #   'diff' = abs(r_llk - mat_llk),
  #   'sig' = abs(r_llk - mat_llk) > 1e-10
  # )
  # 
  # analysis = analysis %>% filter(sig == TRUE)
  # 
  # print(seed)
  # print(analysis)
}

# liklweitz1 parts
consumerData[(consumer_num + 1):(consumer_num + uniCons), 2] = llk
consumer_num = consumer_num + uniCons

# Sum over consumers
llk = -sum(log(consumerData[,2]))

# Check for errors or save output
if (is.na(llk) || llk == Inf || llk == -Inf || !is.numeric(llk) || is.complex(llk)) {
  loglik = 1e+300
} else {
  loglik = llk
}

print(format(loglik, digits=20))

