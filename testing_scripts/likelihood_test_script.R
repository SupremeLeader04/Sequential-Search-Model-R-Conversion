seed = 1
param = c(1,0.7,0.5,0.3,-3)
D = 100
data = read.csv(sprintf("genWeitzDataS%d.csv", seed), header = FALSE)

# Summary statistics storage
summary_stats = list()
summary_stats$b_iH_values = c()
summary_stats$mu_i_H_values = c()
summary_stats$b_ih_values = c()
summary_stats$mu_i_h_values = c()
summary_stats$b_eps_ih_values = c()
summary_stats$eps_i_values = c()
summary_stats$l_stop_values = c()
summary_stats$l_cont_values = c()
summary_stats$l_selection_values = c()
summary_stats$l_choice_values = c()
summary_stats$likelihood_values = c()

# Read random numbers
random_numbers = read.csv('random_numbers_seed1.csv', header = FALSE)
colnames(random_numbers) = c("con", "prod_idx", "draw_idx", "type", "value")

# Extract consumer IDs
consumer = data[,1]
N_cons = length(unique(consumer))

# Extract N_prod (third-to-last column, equivalent to MATLAB's data(:,end-2))
N_prod = data[, ncol(data)-2]
Js = unique(N_prod)
Num_J = length(Js)
consumerData = matrix(0, nrow = N_cons, ncol = 2)
consumer_num = 0
counter = 0
mu_i_sum = 0 # Initialize sum of mu_i means

for (i in 1:Num_J) {
  nalt = Js[i]
  dat = data[N_prod == nalt, , drop = FALSE]
  N_obs = nrow(dat)
  uniCons = N_obs / nalt
  consid2 = matrix(dat[,1], nrow = nalt, ncol = uniCons, byrow = FALSE)
  
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
  for(j in 1:N_cons) {
    last_searched[searched_amt[j], j] = 1
  }
  
  c = exp(param[length(param)]) * rep(1, N_obs)
  X = dat[,4:(3+length(param[1:(length(param)-1)]))]
  xb = rowSums(X * param[1:(length(param)-1)])
  
  # Calculate m
  table = read.csv('tableZ.csv', header=FALSE)
  m = rep(0, N_obs)
  for(j in 1:N_obs) {
    lookupvalue = abs(table[,2] - c[j])
    if(table[1,2] >= c[j] && c[j] >= table[nrow(table),2]) {
      index_m = which.min(lookupvalue)
      m[j] = table[index_m,1]
    } else if(table[1,2] < c[j]) {
      m[j] = -c[j]
    } else if(c[j] < table[nrow(table),2]) {
      m[j] = 4.001
    }
  }
  
  # Load etaDraw from file
  etaDraw = matrix(0, nrow=N_cons, ncol=D)
  eta_subset = random_numbers[random_numbers$type == 1, ]
  for(con in 1:N_cons) {
    con_subset = eta_subset[eta_subset$con == con, ]
    etaDraw[con, ] = con_subset$value[order(con_subset$draw_idx)]
  }
  
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
    
    # Case 2 (Appendix A): consumer searched products and purchased the product searched last
    if(tr == H && H > 1) {
      counter = counter + 1
      
      # Step 1-2: Sampling mu_il
      mu_i_S_bar = matrix(0, nrow=length(S_bar), ncol=D)
      mu_subset = random_numbers[random_numbers$con == con & random_numbers$type == 2, ]
      for(s in 1:length(S_bar)) {
        s_subset = mu_subset[mu_subset$prod_idx == S_bar[s], ]
        mu_i_S_bar[s, ] = s_subset$value[order(s_subset$draw_idx)]
      }
      mu_i[S_bar, ] = mu_i_S_bar
      
      # Step 2: Sampling mu_iH
      upper_bound = rbind(mu_i[1, ] - xb_con_prod[H, con] - mc[H, con],
                          -xb_con_prod[H, con] + xb_con_prod[S_bar, con] + mu_i[S_bar, ])
      b_iH = apply(upper_bound, 2, max)
      uniform_draws_mu_iH = random_numbers[random_numbers$con == con & random_numbers$type == 3 & random_numbers$prod_idx == H, ]$value[order(random_numbers[random_numbers$con == con & random_numbers$type == 3 & random_numbers$prod_idx == H, ]$draw_idx)]
      mu_i[H, ] = qnorm(uniform_draws_mu_iH * (1 - pnorm(b_iH)) + pnorm(b_iH))
      
      # Store random numbers and intermediate values
      summary_stats$b_iH_values = c(summary_stats$b_iH_values, mean(b_iH))
      summary_stats$mu_i_H_values = c(summary_stats$mu_i_H_values, mean(mu_i[H, ]))
      
      # Step 3: Sampling mu_i(1:(H - 1))
      for(h in (H-1):2) {
        b_ih = xb_con_prod[h + 1, con] + mu_i[h + 1, ] - xb_con_prod[h, con]
        uniform_draws_mu_ih = random_numbers[random_numbers$con == con & random_numbers$type == 4 & random_numbers$prod_idx == h, ]$value[order(random_numbers[random_numbers$con == con & random_numbers$type == 4 & random_numbers$prod_idx == h, ]$draw_idx)]
        mu_i[h, ] = qnorm(uniform_draws_mu_ih * (1 - pnorm(b_ih)) + pnorm(b_ih))
        
        # Store random numbers and intermediate values
        summary_stats$b_ih_values = c(summary_stats$b_ih_values, mean(b_ih))
        summary_stats$mu_i_h_values = c(summary_stats$mu_i_h_values, mean(mu_i[h, ]))
      }
      
      # Step 4: Sampling epsilon_i(1:H - 1)
      for(h in 2:(H-1)) {
        b_eps_ih = xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[h, con] - mu_i[h, ]
        uniform_draws_eps_ih = random_numbers[random_numbers$con == con & random_numbers$type == 5 & random_numbers$prod_idx == h, ]$value[order(random_numbers[random_numbers$con == con & random_numbers$type == 5 & random_numbers$prod_idx == h, ]$draw_idx)]
        eps_i[h, ] = qnorm(uniform_draws_eps_ih * pnorm(b_eps_ih))
        
        # Store random numbers and intermediate values
        summary_stats$b_eps_ih_values = c(summary_stats$b_eps_ih_values, mean(b_eps_ih))
        summary_stats$eps_i_values = c(summary_stats$eps_i_values, mean(eps_i[h, ]))
      }
      
      # Stopping rule
      max_values = rbind(mu_i[1, ],
                         xb_con_prod[S_bar, con] + mu_i[S_bar, ] + mc[S_bar, con],
                         mu_i[2:(H-1), ] + xb_con_prod[2:(H-1), con] + eps_i[2:(H-1), ])
      l_stop = 1 - pnorm(apply(max_values, 2, max) - xb_con_prod[H, con] - mu_i[H, ])
      
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
      for(h in 2:(H-1)) {
        l_choice = l_choice * pnorm(xb_con_prod[H, con] + mu_i[H, ] + mc[H, con] - xb_con_prod[h, con] - mu_i[h, ])
      }
      
      # Store rules and likelihood
      summary_stats$l_stop_values = c(summary_stats$l_stop_values,
                                      
                                      mean(l_stop))
      summary_stats$l_cont_values = c(summary_stats$l_cont_values, mean(l_cont))
      summary_stats$l_selection_values = c(summary_stats$l_selection_values, mean(l_selection))
      summary_stats$l_choice_values = c(summary_stats$l_choice_values, mean(l_choice))
      summary_stats$likelihood_values = c(summary_stats$likelihood_values, mean(l_stop * l_cont * l_selection * l_choice))
      
      # Compute mean of mu_i (non-zero elements)
      non_zero_indices = c(1, S_bar, 2:H)  # Indices of mu_i that are populated
      mu_i_non_zero = mu_i[non_zero_indices, ]
      mu_i_mean = mean(mu_i_non_zero)
      mu_i_sum = mu_i_sum + mu_i_mean # Add to running sum
      
    } else {
      l_stop = 1
      l_cont = 1
      l_selection = 1
      l_choice = 1
    }
    
    L_ALL[con, ] = l_stop * l_cont * l_selection * l_choice
  }
  
  llk = rowMeans(L_ALL)
  consumerData[(consumer_num + 1):(consumer_num + uniCons), 1] = t(consid2[1, ])
  consumerData[(consumer_num + 1):(consumer_num + uniCons), 2] = llk
  consumer_num = consumer_num + uniCons
}

# Compute and display overall mean
if(counter > 0) {
  mu_i_overall_mean = mu_i_sum / counter
  cat(sprintf("Overall mean of mu_i for Case 2: %.6f\n", mu_i_overall_mean))
} else {
  cat("No Case 2 consumers found.\n")
}

# Display summary statistics
cat("\n=== R SUMMARY STATISTICS ===\n")
cat(sprintf("Number of Case 2 consumers: %d\n", counter))
cat(sprintf("Mean b_iH: %.6f (std: %.6f)\n", mean(summary_stats$b_iH_values), sd(summary_stats$b_iH_values)))
cat(sprintf("Mean mu_i_H: %.6f (std: %.6f)\n", mean(summary_stats$mu_i_H_values), sd(summary_stats$mu_i_H_values)))
cat(sprintf("Mean b_ih: %.6f (std: %.6f)\n", mean(summary_stats$b_ih_values), sd(summary_stats$b_ih_values)))
cat(sprintf("Mean mu_i_h: %.6f (std: %.6f)\n", mean(summary_stats$mu_i_h_values), sd(summary_stats$mu_i_h_values)))
cat(sprintf("Mean b_eps_ih: %.6f (std: %.6f)\n", mean(summary_stats$b_eps_ih_values), sd(summary_stats$b_eps_ih_values)))
cat(sprintf("Mean eps_i: %.6f (std: %.6f)\n", mean(summary_stats$eps_i_values), sd(summary_stats$eps_i_values)))
cat(sprintf("Mean l_stop: %.6f (std: %.6f)\n", mean(summary_stats$l_stop_values), sd(summary_stats$l_stop_values)))
cat(sprintf("Mean l_cont: %.6f (std: %.6f)\n", mean(summary_stats$l_cont_values), sd(summary_stats$l_cont_values)))
cat(sprintf("Mean l_selection: %.6f (std: %.6f)\n", mean(summary_stats$l_selection_values), sd(summary_stats$l_selection_values)))
cat(sprintf("Mean l_choice: %.6f (std: %.6f)\n", mean(summary_stats$l_choice_values), sd(summary_stats$l_choice_values)))
cat(sprintf("Mean likelihood: %.6f (std: %.6f)\n", mean(summary_stats$likelihood_values), sd(summary_stats$likelihood_values)))

# Sum over consumers
llk = -sum(log(consumerData[,2]))
cat(sprintf("Final log-likelihood: %.10f\n", llk))
cat(sprintf("Counter: %d\n", counter))