# ESTIMATION CODE - GHK - homogeneous parameters - Weitzman, UrsuSeilerHonka 2024

estimate = function(D, seed) {
 
  gc()
  start_time = Sys.time()
  
  ################################################################################
  
  # OPTIONS FOR ESTIMATION #
   
  ################################################################################
  
  optim_options = list(maxit = 6000000, reltol = 1e-6, abstol = 1e-6, fnscale = 1)
  
  ################################################################################
  
  # SETTING PARAMETERS (NOT OBTAINED FROM FILE NAME)
  
  ################################################################################
  
  D = D
  seed = seed
  
  ###############################################################################
  
  # SIMULATION #
  
  ###############################################################################
  
  # simulation inputs
  N_cons = 1000 # num of consumers
  N_prod = 5 # num of products
  param = c(1, 0.7, 0.5, 0.3, -3) # true parameter vector 
                              # [4 brandFE, search cost constant (exp)]
  
  
  # # simulate data
  # source("simWeitz.R")
  # simWeitz(N_cons, N_prod, param, seed)
  
  
  # load simulated data
  data = read.csv(sprintf("genWeitzDataS%d.csv", seed), header = FALSE)

  ################################################################################
  
  # ESTIMATION #
  
  ################################################################################
  
  # initial parameter vector
  param0  = rep(0, length(param))
  
  # do estimation
  source("liklWeitz_ghk_1.R")
  f = function(x) liklWeitz_ghk_1(x, data, D, seed)
  
  # set args to [] as no linear constraints
  A = c()
  b = c()
  Aeq = c()
  beq = c()
  
  # upper and lower bounds on the params
  source("fminsearchcon.R")
  lb = rep(-4, length(param))
  ub = rep(4, length(param))
  nonlcon = c()
  result = fminsearchcon(f, param0, lb, ub, A, b, nonlcon,
                                               optim_options)
  
  be = result$x 
  val = result$fval
  exitflag = result$exitflag 
  output = result$output
  
  # compute standard errors
  source("cal_se_num_hessian.R")
  se = cal_se_num_hessian(be, 0.01, data, D, seed)
  
  # save results
  AS = matrix(c(be, se, val, exitflag), ncol = 1)  # 12x1
  write.csv(AS, sprintf("rezSimWeitz_ghk_D%dS%d.csv", D, seed), row.names = FALSE)
  
  ################################################################################
  
  end_time = Sys.time()
  elapsed_time = end_time - start_time
  
}