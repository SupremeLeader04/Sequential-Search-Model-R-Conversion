liklWeitz_ghk_2 = function(param, dat, D, nalt, seed) {
  
  # Data features
  consumer = dat[,1]
  N_obs = length(consumer)
  N_cons = length(unique(consumer))
  N_prod = nalt
  
  # Choices
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
  
  # Calculate m
  # 1. look-up table method
  table = read.csv('tableZ.csv', header=FALSE)
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
  
  # Set random seed
  set.seed(seed)
  etaDraw = matrix(rnorm(N_cons * D), nrow=N_cons, ncol=D)
  xb_con_prod = matrix(xb, nrow=N_prod, ncol=N_cons)
  mc = matrix(m, nrow=N_prod, ncol=N_cons)
  
  L_ALL = matrix(0, nrow=N_cons, ncol=D)  # initialize LL
  
  for(con in 1:N_cons) {
    H = searched_amt[con]                               # The last searched product
    mu_i = matrix(0, nrow=N_prod, ncol=D)               # Presearch shock vector
    eps_i = matrix(0, nrow=N_prod, ncol=D)              # Postsearch shock vector
    S_bar = setdiff(1:N_prod, 1:H)                      # The set of unsearched products
    tr = tran_ranking[con]                              # Position of purchased product
    
    # Step 1: Sampling mu_i0
    mu_i[1, ] = etaDraw[con, ]  # Sampling for the outside option
    
    # Case 1 (main paper): consumer searched products and only purchased the outside option
    if(tr == 1 && H > 1) {
      
      # Step 2: Sampling mu_iH
      b_iH = mu_i[1, ] - xb_con_prod[H, con] - mc[H, con]
      set.seed(con * seed * 1122)
      mu_i[H, ] = qnorm(runif(length(b_iH)) * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      if (H >= 3) {
        for(h in (H-1):2) {
          b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
          set.seed(con * seed * 2233 + h)
          mu_i[h, ] = qnorm(runif(length(b_ih)) * (1 - pnorm(b_ih)) + pnorm(b_ih))
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
      
      # Step 1-2: Sampling mu_il
      set.seed(seed * con)
      mu_i[S_bar, ] = matrix(rnorm(length(S_bar) * D), nrow=length(S_bar), ncol=D)
      
      # Step 2: Sampling mu_iH
      upper_bound = rbind(mu_i[1, ] - xb_con_prod[H, con] - mc[H, con],
                          -xb_con_prod[H, con] + xb_con_prod[S_bar, con] + mu_i[S_bar, ])
      b_iH = apply(upper_bound, 2, max)
      set.seed(seed * con * 1122)
      mu_i[H, ] = qnorm(runif(length(b_iH)) * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      if (H >= 3) {
        for(h in (H-1):2) {
          b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
          set.seed(con * seed * 2233 + h)
          mu_i[h, ] = qnorm(runif(length(b_ih)) * (1 - pnorm(b_ih)) + pnorm(b_ih))
        }
      }
      
      # Step 4: Sampling epsilon_i(1:H - 1)
      if (H >= 3) {
        for(h in 2:(H-1)) {
          b_eps_ih = xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[h, con] - mu_i[h, ]
          set.seed(seed * con * 3344 + h)
          eps_i[h, ] = qnorm(runif(length(b_eps_ih)) * pnorm(b_eps_ih))
        }
      }
      
      # Stopping rule
      if (H >= 3) {
        middle_term = mu_i[2:(H-1), ] + xb_con_prod[2:(H-1), con] + eps_i[2:(H-1), ]
        l_stop <- 1 - pnorm(apply(rbind(mu_i[1, ], 
                                        (xb_con_prod[S_bar, con] + mu_i[S_bar, ] + mc[S_bar, con]), 
                                        middle_term), 
                                  2, max) - xb_con_prod[H, con] - mu_i[H, ])
      } else {
        l_stop = 1 - pnorm(apply(rbind(mu_i[1, ], 
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
      for(h in 2:(H-1)) {
        l_choice = l_choice * pnorm(xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[h, con] - mu_i[h, ])
      }
      
      # Case 3 (Appendix A): consumer searched products and purchased a product neither not searched last nor outside option
    } else if((tr > 1) && (tr < H)) {
      
      # Step 1-2: Sampling mu_il
      set.seed(seed * con)
      mu_i[S_bar, ] = matrix(rnorm(length(S_bar) * D), nrow=length(S_bar), ncol=D)
      
      # Step 2: Sampling mu_iH
      upper_bound = rbind(mu_i[1, ] - xb_con_prod[H, con] - mc[H, con],
                          -xb_con_prod[H, con] + xb_con_prod[S_bar, con] + mu_i[S_bar, ])
      b_iH = apply(upper_bound, 2, max)
      set.seed(seed * con * 1122)
      mu_i[H, ] = qnorm(runif(length(b_iH)) * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      for(h in (H-1):2) {
        b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
        set.seed(con * seed * 2233 + h)
        mu_i[h, ] = qnorm(runif(length(b_ih)) * (1 - pnorm(b_ih)) + pnorm(b_ih))
      }
      
      # Step 4: Sampling epsilon for the purchased option
      b_eps_upper = xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[tr, con] - mu_i[tr, ]
      max_values = rbind(xb_con_prod[S_bar, con] + mu_i[S_bar, ] + mc[S_bar, con],
                         mu_i[1, ])
      b_eps_lower = apply(max_values, 2, max) - xb_con_prod[tr, con] - mu_i[tr, ]
      set.seed(seed * con * 3344)
      eps_i[tr, ] = qnorm(runif(length(b_eps_upper)) * (pnorm(b_eps_upper) - pnorm(b_eps_lower)) + pnorm(b_eps_lower))
      
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
  return(llk)
}