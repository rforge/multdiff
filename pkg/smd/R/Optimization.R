## Gauss-Newton type algorithm:
## 'f' computes the least squares loss function
## 'gr' computes the gradient of the least squares loss
## 'quad' computes, as a function of an index along the beta vector and
## beta, the quadratic coefficient for the quadratic approximation to
## the loss function in the corresponding direction.
## 'lambda' is a vector of penalty parameters, which will be used in
## decreasing order, and for each new lambda the parameters estimated
## for the previous lambda will be used as a warm start.
## The function returns a list with the vector of lambda values as the first entry and the estimated beta parameters as a matrix in the second entry. Each column in the matrix corresponds to one lambda value.

GNcoordinateDescent <- function(beta, f, gr, quad,
                                lambda = 1,
                                reltol = sqrt(.Machine$double.eps),
                                trace = 0) {
  if(trace > 0)
    cat("lambda\tdf\tpenalty\t\tloss+penalty\n")
  beta <- matrix(beta, nrow = length(beta), ncol = length(lambda))
  pIndex <- seq_len(dim(beta)[1])
  lambda <- sort(lambda, decreasing = TRUE)
  for(j in seq_along(lambda)){
    if(j > 1) 
      beta[ , j] <- beta[ , j - 1]
    pen <- lambda[j]*sum(abs(beta[ , j]))
    val <-  f(beta[ , j]) + pen
    grad = gr(beta[ , j])
    repeat {
      oldval <- val
      for(i in pIndex) {
        if(beta[i, j] == 0) {
          if(abs(grad[i]) > lambda[j]) {
            beta0 <- beta[ , j]
            alpha <-  quad(i, beta[ , j])
            if(grad[i] < -lambda[j]) {
              beta[i, j] <- -(lambda[j] + grad[i])/(2*alpha)
            } else {  ## grad[i] > lambda[j]
              beta[i, j] <- (lambda[j] - grad[i])/(2*alpha)
            }
            pen <- pen + lambda[j]*abs(beta[i, j])
            grad <- gr(beta[ , j])
          }
        } else { ## beta[i, j] != 0
          beta0 <- beta[ , j]
          beta0i <- beta[i , j]
          alpha <-  quad(i, beta[ , j])
          di <- -2*alpha*beta0i + grad[i]
          if(abs(di) <= lambda[j]) {
            beta[i, j] <- 0
          } else { ## abs(di) > lambda[j]
            if(di < -lambda[j]) {
              beta[i, j] <- beta0i - (lambda[j] + grad[i])/(2*alpha)
            } else { ## di > lambda[j]
              beta[i, j] <- beta0i + (lambda[j] - grad[i])/(2*alpha)
            }
          }
          pen <- pen + lambda[j]*(abs(beta[i, j]) - abs(beta0i))
          grad <- gr(beta[ , j])
        }
      }
      val <- f(beta[ , j]) + pen
      if(oldval-val < reltol*(abs(val) + reltol))
        break
      if(trace == 2)
        cat(lambda[j], "\t", sum(abs(beta[,j]) > 0), "\t", val, "\n")
    }
    if(trace == 1)
      cat(lambda[j], "\t", sum(abs(beta[,j]) > 0), "\t", pen, "\t", val, "\n")
  }
  return(list(lambda = lambda, beta = beta))
}

## Non-optimal solution based on 'optimize'. This solution solves the
## coodinate wise optimizations in a generic way based on the
## 'optimize' function. It is relatively slow.

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
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) + lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(0, gamma0))$minimum
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) - lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(-gamma0, 0))$minimum
            }
            pen <- pen + lambda[j]*abs(beta[i, j])
            val <- f(beta[ , j]) + pen
          }
        } else {
          beta0 <- beta[ , j]
          beta0i <- beta[i , j]
          if(beta0i > 0) {
             tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) + lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(0, beta0i+gamma0))$minimum
           } else { ## beta0i < 0
            tmpf <- function(b) {
              beta0[i] <- b
              f(beta0) - lambda[j]*b
            }
            beta[i, j] <- optimize(tmpf, c(beta0i-gamma0, 0))$minimum
          }
          pen <- pen + lambda[j]*(abs(beta[i, j]) - abs(beta0i))
          val <- f(beta[ , j]) + pen
          grad <- gr(beta[ , j])
        }
      }
      if(oldval-val < reltol*(abs(val) + reltol))
        break
      cat(j, lambda[j], sum(abs(beta[,j]) > 0), f(beta[ ,j]) + lambda[j]*sum(abs(beta)), val, "\n")
    }
  }
  return(list(lambda = lambda, beta = beta))
}

## WARNING: This coordinate descent algorithm was experimental,
## and in its current form not stable! It relies on coordinate wise
## computations of quadratic approximations and then minimization of
## the quadratic function, but the steps used to make the quadratic
## approximation are not adaptive, which do make the algorithm diverge
## for some cases. The algorithm is general, and requires only the
## function and gradient. The other algorithm implemented above requires
## a separate function to compute the quadratic approximation, which is
## well suited for non-linear least squares problems. 

ObsoleteCoordinateDescent <- function(beta, f, gr, lambda = 1, reltol = sqrt(.Machine$double.eps), epsilon = 0, gamma = 1e-2) {
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
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) + lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(0, gamma0))$minimum
              ## Dubious quadratic approximation step
              ## beta0[i] <- gamma0
              ## newval <- f(beta0) + pen 
              ## alpha <-  2*((newval - val)/gamma0^2 - grad[i]/gamma0)
              ## beta[i, j] <- -(lambda[j] + grad[i])/alpha
            } else {  ## grad[i] > lambda[j]
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) - lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(-gamma0, 0))$minimum
              ## Dubious quadratic approximation step
              ## beta0[i] <- - gamma0
              ## newval <- f(beta0) + pen
              ## alpha <-  2*((newval - val)/gamma0^2 + grad[i]/gamma0)
              ## beta[i, j] <- (lambda[j] - grad[i])/alpha
            }
            pen <- pen + lambda[j]*abs(beta[i, j])
            val <- f(beta[ , j]) + pen
          }
        } else {
          beta0 <- beta[ , j]
          beta0i <- beta[i , j]
          if(beta0i > 0) {
             tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) + lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(0, beta0i+gamma0))$minimum
             ## Dubious quadratic approximation step
             ## if(grad[i] < -lambda[j]) {
             ## beta0[i] <- beta0i + gamma0
             ## newval <- f(beta0) + pen 
             ## alpha <-  2*((newval - val)/gamma0^2 - grad[i]/gamma0)
             ## beta[i, j] <- beta[i ,j] - (lambda[j] + grad[i])/alpha
             ## pen <- pen + lambda[j]*(beta[i, j] - beta0i)
             ## val <- f(beta[ , j]) + pen 
             ## } else { ## grad[i] > -lambda[j]
             ## beta0[i] <- 0
             ## newval <- f(beta0) + pen 
             ## if(newval < val + lambda[j]*beta[i, j]) {
             ##   pen <- pen - lambda[j]*beta[i, j]
             ##   beta[i, j] <- 0
             ##   val <- newval - lambda[j]*beta[i, j]
             ## } else {
             ##   if(beta[i, j] < gamma0) {
             ##     alpha <-  2*((newval - val)/beta[i, j]^2 + grad[i]/beta[i, j])
             ##   } else {
             ##     beta0[i] <- beta[i, j] - gamma0
             ##     newval <- f(beta0) + pen 
             ##     alpha <-  2*((newval - val)/gamma0^2 - grad[i]/gamma0)
             ##   }
             ##   beta[i, j] <- beta[i, j] - (lambda[j] + grad[i])/alpha
             ##   pen <- pen + lambda[j]*(beta[i, j] - beta0i)
             ##   val <- f(beta[ , j]) + pen
             ## }
             ## }
          } else { ## beta0i < 0
            tmpf <- function(b) {
              beta0[i] <- b
              f(beta0) - lambda[j]*b
            }
            beta[i, j] <- optimize(tmpf, c(beta0i-gamma0, 0))$minimum
            ## if(grad[i] > lambda[j]) {
            ##   beta0[i] <- beta0i - gamma0
            ##   newval <- f(beta0) + pen 
            ##   alpha <-  2*((newval - val)/gamma0^2 + grad[i]/gamma0)
            ##   beta[i, j] <- beta[i, j] + (lambda[j] - grad[i])/alpha
            ##   pen <- pen + lambda[j]*(beta0i - beta[i ,j])
            ##   val <- f(beta[ , j]) + pen 
            ## } else { ## grad[i] < lambda[j]
            ##   beta0[i] <- 0
            ##   newval <- f(beta0) + pen 
            ##   if(newval < val - lambda[j]*beta[i ,j]) {
            ##     pen <- pen + lambda[j]*beta[i, j]
            ##     beta[i, j] <- 0
            ##     val <- newval + lambda[j]*beta[i ,j]
            ##   } else {
            ##     if(beta[i ,j] > -gamma0) {
            ##       alpha <-  2*((newval - val)/beta[i ,j]^2 - grad[i]/beta[i ,j])
            ##     } else {
            ##       beta0[i] <- beta[i ,j] + gamma0
            ##       newval <- f(beta0) + pen
            ##       alpha <-  2*((newval - val)/gamma0^2 + grad[i]/gamma0)
            ##     }
            ##     beta[i ,j] <- beta[i ,j] + (lambda[j] - grad[i])/alpha
            ##     pen <- pen + lambda[j]*(beta0i - beta[i, j])
            ##     newval <- f(beta[ , j]) + pen
            ##   }
            ##  }
          }
          pen <- pen + lambda[j]*(abs(beta[i, j]) - abs(beta0i))
          val <- f(beta[ , j]) + pen
          grad <- gr(beta[ , j])
        }
      }
      if(oldval-val < reltol*(abs(val) + reltol))
        break
      cat(j, lambda[j], sum(abs(beta[,j]) > 0), f(beta[ ,j]) + lambda[j]*sum(abs(beta)), val, "\n")
    }
  }
  return(list(lambda = lambda, beta = beta))
}
