## WARNING: This coordinate descent algorithm is highly experimental,
## and in its current form not stable! It relies on coordinate wise
## computations of quadratic approximations and then minimization of
## the quadratic function, but the steps used to make the quadratic
## approximation are not adaptive, which do make the algorithm diverge
## for some cases. The algorithm is general, and requires only the
## function and gradient. Another algorithm will be implemented that
## is specific to the non-linear least squares problem.

coordinateDescent <- function(beta, f, gr, lambda = 1, reltol = sqrt(.Machine$double.eps), epsilon = 0, gamma = 1e-2) {
  gamma0 <- gamma
  beta <- matrix(beta, nrow = length(beta), ncol = length(lambda))
  pIndex <- seq_len(dim(beta)[1])
  for(j in seq_along(lambda)){
    pen <- lambda[j]*sum(abs(beta[ , j]))
    if(j > 1) 
      beta[ , j] <- beta[ , j - 1]
    val <-  f(beta[ , j]) + pen
    grad = gr(beta[ , j])
    repeat {
      oldval <- val
      for(i in pIndex) {
        if(beta[i, j] == 0) {
          if(abs(grad[i]) > lambda[j] + epsilon) {
            beta0 <- beta[ , j]
            if(grad[i] < -lambda[j]) {
              beta0[i] <- gamma0
              newval <- f(beta0) + pen 
              alpha <-  2*((newval - val)/gamma0^2 - grad[i]/gamma0)
              beta[i, j] <- -(lambda[j] + grad[i])/alpha
            } else {  ## grad[i] > lambda[j]
              beta0[i] <- - gamma0
              newval <- f(beta0) + pen
              alpha <-  2*((newval - val)/gamma0^2 + grad[i]/gamma0)
              beta[i, j] <- (lambda[j] - grad[i])/alpha
            }
            pen <- pen + lambda[j]*abs(beta[i, j])
            val <- f(beta[ , j]) + pen
          }
        } else {
          beta0 <- beta[ , j]
          beta0i <- beta[i , j]
          if(beta0i > 0) {
            if(grad[i] < -lambda[j]) {
              beta0[i] <- beta0i + gamma0
              newval <- f(beta0) + pen 
              alpha <-  2*((newval - val)/gamma0^2 - grad[i]/gamma0)
              beta[i, j] <- beta[i ,j] - (lambda[j] + grad[i])/alpha
              pen <- pen + lambda[j]*(beta[i, j] - beta0i)
              val <- f(beta[ , j]) + pen 
            } else { ## grad[i] > -lambda[j]
              beta0[i] <- 0
              newval <- f(beta0) + pen 
              if(newval < val + lambda[j]*beta[i, j]) {
                pen <- pen - lambda[j]*beta[i, j]
                beta[i, j] <- 0
                val <- newval - lambda[j]*beta[i, j]
              } else {
                if(beta[i, j] < gamma0) {
                  alpha <-  2*((newval - val)/beta[i, j]^2 + grad[i]/beta[i, j])
                } else {
                  beta0[i] <- beta[i, j] - gamma0
                  newval <- f(beta0) + pen 
                  alpha <-  2*((newval - val)/gamma0^2 - grad[i]/gamma0)
                }
                beta[i, j] <- beta[i, j] - (lambda[j] + grad[i])/alpha
                pen <- pen + lambda[j]*(beta[i, j] - beta0i)
                val <- f(beta[ , j]) + pen
              }
            }
          } else { ## beta0i < 0
            if(grad[i] > lambda[j]) {
              beta0[i] <- beta0i - gamma0
              newval <- f(beta0) + pen 
              alpha <-  2*((newval - val)/gamma0^2 + grad[i]/gamma0)
              beta[i, j] <- beta[i, j] + (lambda[j] - grad[i])/alpha
              pen <- pen + lambda[j]*(beta0i - beta[i ,j])
              val <- f(beta[ , j]) + pen 
            } else { ## grad[i] < lambda[j]
              beta0[i] <- 0
              newval <- f(beta0) + pen 
              if(newval < val - lambda[j]*beta[i ,j]) {
                pen <- pen + lambda[j]*beta[i, j]
                beta[i, j] <- 0
                val <- newval + lambda[j]*beta[i ,j]
              } else {
                if(beta[i ,j] > -gamma0) {
                  alpha <-  2*((newval - val)/beta[i ,j]^2 - grad[i]/beta[i ,j])
                } else {
                  beta0[i] <- beta[i ,j] + gamma0
                  newval <- f(beta0) + pen
                  alpha <-  2*((newval - val)/gamma0^2 + grad[i]/gamma0)
                }
                beta[i ,j] <- beta[i ,j] + (lambda[j] - grad[i])/alpha
                pen <- pen + lambda[j]*(beta0i - beta[i, j])
                newval <- f(beta[ , j]) + pen
              }
            }
          }
          grad <- gr(beta[ , j])
        }
      }
      if(abs(oldval-val) < reltol*(abs(val) + reltol))
        break
      cat(j, lambda[j], sum(abs(beta[,j]) > 0), f(beta[ ,j]) + lambda[j]*sum(abs(beta)), val, "\n")
    }
  }
  return(list(lambda = lambda, beta = beta))
}
