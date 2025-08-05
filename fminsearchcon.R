fminsearchcon <- function(fun, x0, LB = NULL, UB = NULL, A = NULL, b = NULL, 
                          nonlcon = NULL, options = NULL, ...) {
  # FMINSEARCHCON: Extension of Nelder-Mead with general inequality constraints
  # R port of the MATLAB function by John D'Errico
  
  # Store original size of x0
  xsize <- dim(x0)
  if (is.null(xsize)) xsize <- length(x0)
  x0 <- as.vector(x0)
  n <- length(x0)
  
  # Handle bounds
  if (is.null(LB) || length(LB) == 0) {
    LB <- rep(-Inf, n)
  } else {
    LB <- as.vector(LB)
  }
  
  if (is.null(UB) || length(UB) == 0) {
    UB <- rep(Inf, n)
  } else {
    UB <- as.vector(UB)
  }
  
  # Size checks
  if (n != length(LB) || n != length(UB)) {
    stop("x0 is incompatible in size with either LB or UB.")
  }
  
  # Handle A and b
  if (is.null(A)) A <- matrix(nrow = 0, ncol = n)
  if (is.null(b)) b <- numeric(0)
  
  if ((nrow(A) == 0 && length(b) > 0) || (length(b) == 0 && nrow(A) > 0)) {
    stop("Sizes of A and b are incompatible")
  }
  
  if (nrow(A) > 0) {
    A <- as.matrix(A)
    b <- as.vector(b)
    if (nrow(A) != length(b)) {
      stop("Sizes of A and b are incompatible")
    }
    if (ncol(A) != n) {
      stop("A is incompatible in size with x0")
    }
  }
  
  # Test feasibility of initial values
  if (nrow(A) > 0) {
    if (any(A %*% x0 > b)) {
      stop("Infeasible starting values (linear inequalities failed).")
    }
  }
  
  if (!is.null(nonlcon)) {
    x0_reshaped <- if (is.null(dim(xsize))) x0 else array(x0, dim = xsize)
    cons <- nonlcon(x0_reshaped, ...)
    if (any(cons > 0)) {
      stop("Infeasible starting values (nonlinear inequalities failed).")
    }
  }
  
  # Set default options
  if (is.null(options)) {
    options <- list(
      reltol = 1e-4,
      abstol = 1e-8,
      maxit = 1000,
      trace = 0
    )
  }
  
  # Create parameter structure
  params <- list(
    args = list(...),
    LB = LB,
    UB = UB,
    fun = fun,
    n = n,
    xsize = xsize,
    A = A,
    b = b,
    nonlcon = nonlcon
  )
  
  # Classify bounds
  # 0 --> unconstrained variable
  # 1 --> lower bound only
  # 2 --> upper bound only
  # 3 --> dual finite bounds
  # 4 --> fixed variable
  params$BoundClass <- numeric(n)
  for (i in 1:n) {
    k <- as.numeric(is.finite(LB[i])) + 2 * as.numeric(is.finite(UB[i]))
    params$BoundClass[i] <- k
    if (k == 3 && LB[i] == UB[i]) {
      params$BoundClass[i] <- 4
    }
  }
  
  # Transform starting values
  x0u <- numeric(0)
  k <- 1
  for (i in 1:n) {
    switch(as.character(params$BoundClass[i]),
           "1" = {
             # lower bound only
             if (x0[i] <= LB[i]) {
               x0u[k] <- 0
             } else {
               x0u[k] <- sqrt(x0[i] - LB[i])
             }
             k <- k + 1
           },
           "2" = {
             # upper bound only
             if (x0[i] >= UB[i]) {
               x0u[k] <- 0
             } else {
               x0u[k] <- sqrt(UB[i] - x0[i])
             }
             k <- k + 1
           },
           "3" = {
             # lower and upper bounds
             if (x0[i] <= LB[i]) {
               x0u[k] <- -pi/2
             } else if (x0[i] >= UB[i]) {
               x0u[k] <- pi/2
             } else {
               x0u[k] <- 2 * (x0[i] - LB[i]) / (UB[i] - LB[i]) - 1
               x0u[k] <- 2 * pi + asin(max(-1, min(1, x0u[k])))
             }
             k <- k + 1
           },
           "0" = {
             # unconstrained variable
             x0u[k] <- x0[i]
             k <- k + 1
           },
           "4" = {
             # fixed variable - don't increment k
           }
    )
  }
  
  # Check if all variables were fixed
  if (length(x0u) == 0) {
    x <- xtransform(x0u, params)
    x <- if (is.null(dim(xsize))) x else array(x, dim = xsize)
    fval <- do.call(params$fun, c(list(x), params$args))
    
    output <- list(
      iterations = 0,
      funcCount = 1,
      algorithm = "nelder-mead",
      message = "All variables were held fixed by the applied bounds"
    )
    
    return(list(x = x, fval = fval, exitflag = 0, output = output))
  }
  
  # Define the internal objective function
  intrafun <- function(x, params) {
    # Transform variables
    xtrans <- xtransform(x, params)
    
    # Test linear constraints
    if (nrow(params$A) > 0) {
      if (any(params$A %*% xtrans > params$b)) {
        return(Inf)
      }
    }
    
    # Reshape for function calls
    xtrans_shaped <- if (is.null(dim(params$xsize))) xtrans else array(xtrans, dim = params$xsize)
    
    # Test nonlinear constraints
    if (!is.null(params$nonlcon)) {
      cons <- do.call(params$nonlcon, c(list(xtrans_shaped), params$args))
      if (any(cons > 0)) {
        return(Inf)
      }
    }
    
    # Evaluate objective function
    fval <- do.call(params$fun, c(list(xtrans_shaped), params$args))
    return(fval)
  }
  
  # Create control list for optim
  control <- list(
    reltol = ifelse(is.null(options$reltol), 1e-4, options$reltol),
    abstol = ifelse(is.null(options$abstol), 1e-8, options$abstol),
    maxit = ifelse(is.null(options$maxit), 1000, options$maxit),
    trace = ifelse(is.null(options$trace), 0, options$trace)
  )
  
  # Call optim with Nelder-Mead
  result <- optim(
    par = x0u,
    fn = intrafun,
    params = params,
    method = "Nelder-Mead",
    control = control
  )
  
  # Transform back to original space
  x <- xtransform(result$par, params)
  x <- if (is.null(dim(xsize))) x else array(x, dim = xsize)
  
  # Determine exit flag
  exitflag <- switch(as.character(result$convergence),
                     "0" = 1,   # successful convergence
                     "1" = 0,   # iteration limit reached
                     "10" = -1, # degeneracy in Nelder-Mead
                     0          # other
  )
  
  # Create output structure
  output <- list(
    iterations = result$counts[1],
    funcCount = result$counts[2],
    algorithm = "nelder-mead",
    message = result$message
  )
  
  return(list(x = x, fval = result$value, exitflag = exitflag, output = output))
}

# Helper function to transform variables back to original space
xtransform <- function(x, params) {
  xtrans <- numeric(params$n)
  k <- 1
  
  for (i in 1:params$n) {
    switch(as.character(params$BoundClass[i]),
           "1" = {
             # lower bound only
             xtrans[i] <- params$LB[i] + x[k]^2
             k <- k + 1
           },
           "2" = {
             # upper bound only
             xtrans[i] <- params$UB[i] - x[k]^2
             k <- k + 1
           },
           "3" = {
             # lower and upper bounds
             xtrans[i] <- (sin(x[k]) + 1) / 2
             xtrans[i] <- xtrans[i] * (params$UB[i] - params$LB[i]) + params$LB[i]
             # Ensure within bounds due to floating point issues
             xtrans[i] <- max(params$LB[i], min(params$UB[i], xtrans[i]))
             k <- k + 1
           },
           "4" = {
             # fixed variable
             xtrans[i] <- params$LB[i]
           },
           "0" = {
             # unconstrained variable
             xtrans[i] <- x[k]
             k <- k + 1
           }
    )
  }
  
  return(xtrans)
}