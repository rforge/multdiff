### R code for optimizing non-convex l1-penalized loss functions.
###
###     Copyright (C) 2014 Niels Richard Hansen.
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2, or (at your option) any
### later version.
###
### These functions are distributed in the hope that they will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/

## Misc. helper functions. Not exported.

penalty <- function(beta, lambda) 
  sum(lambda * abs(beta))

##' Optimize l1-penalized loss using a coordinate wise descent algorithm
##' 
##' The function computes a matrix of parameter estimates. Each column corresponds to a
##' value of the penalty parameter in \code{lambda}. The estimates are computed in decreasing order 
##' of the penalty parameters, and for each column the previous is used as a warm start.
##' 
##' The algorithm relies on iterative optimization of an l1-penalized 
##' quadratic approximation of the loss using a standard coordinate wise descent algorithm. An 
##' outer backtracking step is added to ensure that the algorithm takes descent steps. 
##' 
##' This function relies on three auxiliary functions. The loss function \code{f}, its gradient \code{gr}
##' and a third function, \code{quad}, that computes the coefficient of the quadratic approximation 
##' for the coordinate wise optimization. 
##' 
##' The function returns a list with the vector of lambda values as the first entry 
##' and the estimated beta parameters as a matrix in the second entry. Each column 
## in the matrix corresponds to one lambda value.
##' 
##' @param beta a \code{numeric}. The vector of initial parameter values.
##' @param f a \code{function}. The loss function. A function of \code{beta}.
##' @param gr a \code{function}. The gradient of the loss. A function of \code{beta}.
##' @param quad a \code{function}. The coordinate wise quadratic approximation term. 
##'         A function of the coordinate index and \code{beta}.
##' @param lambda a \code{numeric}. A sequence of penalties. The default value \code{NULL} 
##'         implies an automatic computation of a suitable sequence.
##' @param penalty.factor a \code{numeric}. A vector of weight factors.
##' @param rho a \code{numeric}. Step length control in the backtracking.
##' @param c a \code{numeric}. Sufficient decrease control in the backtracking.
##' @param reltol a \code{numeric}. Controls the convergence criterion.
##' @param trace a \code{numeric}. Values above 0 prints increasing amounts of trace information.
##' @return A \code{list} of length 2. The first entry contains \code{lambda} and the second the
##'         matrix of parameter estimates. Each column in the matrix corresponds to an entry in \code{lambda}.
##' @export

coordinateDescentQuad <- function(beta, 
                                  f, 
                                  gr, 
                                  quad,
                                  lambda = 1,
                                  penalty.factor = rep(1, p),
                                  rho = 0.9,
                                  c = 1e-4,
                                  reltol = sqrt(.Machine$double.eps),
                                  trace = 0
                                  ) {
  grtol <- reltol ## TODO: Does grtol need to be a different numeric value?
  if(trace > 0)
    message("lambda\tdf\tpenalty\t\tloss+penalty")
  beta <- matrix(beta, nrow = length(beta), ncol = length(lambda))
  pIndex <- seq_len(nrow(beta))
  lambda <- sort(lambda, decreasing = TRUE)
  for(j in seq_along(lambda)){
    if(j > 1) 
      beta[ , j] <- beta[ , j - 1]
    lamb <- lambda[j] * penalty.factor
    pen <- penalty(beta[, j], lamb)
    val <-  f(beta[ , j]) + pen
    grad = gr(beta[ , j])
    repeat {
      oldval <- val
      for(i in pIndex) {
        if(abs(grad[i]) > grtol) {
          if(beta[i, j] == 0) {
            if(abs(grad[i]) > lamb[i]) {
              alpha <-  quad(i, beta[, j])
              if(grad[i] < - lamb[i]) {
                beta[i, j] <- - (lamb[i] + grad[i]) / (2 * alpha)
                cg <- c * (grad[i] + lamb[i])
              } else {  ## grad[i] > lamb[i]
                beta[i, j] <- (lamb[i] - grad[i]) / (2 * alpha)
                cg <- c * (grad[i] - lamb[i])
              }
              ## Backtracking
              val <- f(beta[, j]) + pen + lamb[i] * abs(beta[i, j])
              while(val > oldval + beta[i, j] * cg) {
                beta[i, j] <- rho * beta[i, j]
                val <- f(beta[, j]) + pen + lamb[i] * abs(beta[i, j])
              }
              pen <- penalty(beta[, j], lamb)
              grad <- gr(beta[, j])
            }
          } else { ## beta[i, j] != 0
            beta0i <- beta[i, j]
            pen <- pen - lamb[i] * abs(beta0i)
            cg <- c * (grad[i] + sign(beta0i) * lamb[i])
            alpha <-  quad(i, beta[, j])
            di <- 2 * alpha * beta0i - grad[i]
            if(abs(di) <= lamb[i]) {
              beta[i, j] <- 0
            } else { ## abs(di) > lamb[i]
              if(di < - lamb[i]) {
                delta <- (lamb[i] - grad[i]) / (2 * alpha)
              } else { ## di > lamb
                delta <- - (lamb[i] + grad[i]) / (2 * alpha)
              }
              ## Backtracking
              beta[i, j] <- beta0i + delta
              val <- f(beta[, j]) + pen + lamb[i] * abs(beta[i, j])
              while(val > oldval + delta * cg) {
                delta <- rho * delta
                beta[i, j] <- beta0i + delta
                val <- f(beta[, j]) + pen + lamb[i] * abs(beta[i, j])
              }
            }
            pen <- penalty(beta[, j], lamb)
            grad <- gr(beta[, j])
          }
        }
      }
      val <- f(beta[, j]) + pen
      if(oldval - val < reltol * (abs(val) + reltol))
        break
      if(trace == 2)
        message(signif(lambda[j], 2), 
                "\t", sum(abs(beta[,j]) > 0), 
                "\t", signif(pen, 2), 
                "\t", signif(val, 2))
    }
    if(trace == 1)
      message(signif(lambda[j], 2), 
              "\t", sum(abs(beta[,j]) > 0), 
              "\t", signif(pen, 2), 
              "\t", signif(val, 2))
  }
  return(list(lambda = lambda, beta = beta))
}


## Coordinate descent algorithm for non-linear least squares

coordinateDescent <- function(beta, y, xi, dxi, 
                              rho = 0.9,
                              c = 1e-4,
                              lambda = 1,
                              reltol = sqrt(.Machine$double.eps),
                              trace = 0, 
                              penalty.factor = rep(1, p),
                              N = 1000) {
  kk <- k <- 1
  eps <- .Machine$double.eps
  if(trace > 0)
    message("lambda\tdf\tpenalty\t\tloss+penalty")
  p <- length(beta)
  beta <- matrix(beta, nrow = p, ncol = length(lambda))
  penalty.factor <- penalty.factor * length(y)
  pIndex <- seq_len(p)
  delta <- numeric(p)
  beta0 <- numeric(p)
  deltak <- 1
  lambda <- sort(lambda, decreasing = TRUE)
  for(j in seq_along(lambda)){
    if(j > 1) 
      beta[, j] <- beta[, j - 1]
    lamb <- lambda[j] * penalty.factor
    pen <- penalty(beta[, j], lamb)
    r <- xi(beta[, j]) - y
    ell <-  drop(crossprod(r)) 
    loss <- ell + pen
    repeat {
      oldloss <- loss
      beta0 <- beta[, j]
      delta <- numeric(p)
      X <- dxi(beta0) ## n x p matrix
      cpX <- vector("list", p)
      gamma <- 2 * drop(crossprod(X, r))
      sqloss <- penalty(beta0, lamb) 
      repeat {
        sqoldloss <- sqloss
        for (i in pIndex) {
          if (beta0[i] == 0) {
            if (abs(gamma[i]) > lamb[i]) {
              alpha <- drop(crossprod(X[, i]))
              if (gamma[i] < - lamb[i]) {
                d <- - (lamb[i] + gamma[i]) / (2 * alpha)
              } else {  ## gamma > lambda[j]
                d <- (lamb[i] - gamma[i]) / (2 * alpha)
              }
              if (is.null(cpX[[i]]))
                cpX[[i]] <- drop(crossprod(X, X[, i]))
              gamma <- gamma + (2 * d) * cpX[[i]]
              delta[i] <- delta[i] + d
              beta0[i] <- beta0[i] + d
            }
          } else { ## beta0[i] != 0 
            alpha <- drop(crossprod(X[, i]))
            di <- 2 * alpha * beta0[i] - gamma[i]
            if (abs(di) <= lamb[i]) {
              d <- - beta0[i]
            } else { ## abs(di) > lambda[j]
              if (di < - lamb[i]) {
                d <- (lamb[i] - gamma[i]) / (2 * alpha)
              } else { ## di > lambda[j]
                d <- - (lamb[i] + gamma[i]) / (2 * alpha)
              }
            }
            if (is.null(cpX[[i]]))
              cpX[[i]] <- drop(crossprod(X, X[, i]))
            gamma <- gamma + (2 * d) * cpX[[i]]
            delta[i] <- delta[i] + d
            beta0[i] <- beta0[i] + d
          }
        }
        pen <- penalty(beta0, lamb) 
        Xdelta <- X %*% delta  ## This can be done more efficiently due to the zeroes in delta
        sqloss <- sum(Xdelta^2) + 2 * drop(crossprod(r, Xdelta)) + pen
        if(trace >= 4)
            message("Quadratic loop:", 
                    signif(lambda[j], 2), 
                    "\t", sum(abs(beta0) > 0), 
                    "\t", signif(pen, 2), 
                    "\t", signif(sqloss, 2))
        if(sqoldloss - sqloss < reltol * (abs(sqloss) + reltol))
          break
        if (kk > N) {
          message("Max iterations exit from inner loop.")
          kk <- 1
          break
        }
        kk <- kk + 1
      }
      ## Backtracking
      oldpen <- penalty(beta[, j], lamb)
      cg <- c * (2 * drop(crossprod(r, Xdelta)) + pen - oldpen)
#       deltak <- 1
#       r <- xi(beta0) - y
#       ell <-  drop(crossprod(r))
#       loss <- ell + pen
      beta0 <- beta[, j] + deltak * delta
      r <- xi(beta0) - y
      ell <-  drop(crossprod(r))
      pen <- penalty(beta0, lamb)
      loss <- ell + pen
      while(loss > oldloss + deltak * cg) {
        if(trace >= 3)
          message("Backtrack (inner): ", 
                  signif(lambda[j], 2), 
                  "\t", sum(abs(beta0) > 0), 
                  "\t", signif(pen, 2), 
                  "\t", signif(loss, 2))
        deltak <- rho * deltak
        if (deltak < eps) {
          deltak <- 1
          break
        }
        beta0 <- beta[, j] + deltak * delta
        r <- xi(beta0) - y
        ell <-  drop(crossprod(r))
        pen <- penalty(beta0, lamb)
        loss <- ell + pen
      }
      beta[, j] <- beta0
      if(trace >= 2)
        message("Backtrack (outer): ", 
                signif(lambda[j], 2), 
                "\t", sum(abs(beta[, j]) > 0), 
                "\t", signif(pen, 2), 
                "\t", signif(loss, 2))
      if(oldloss - loss < reltol * (abs(loss) + reltol))
        break
      if (k > N) {
        message("Max iterations exit in outer loop.")
        k <- 1
        break
      }
      k <- k + 1
    }
    if(trace == 1)
      message(signif(lambda[j], 2), 
              "\t", sum(abs(beta[, j]) > 0), 
              "\t", signif(pen, 2), "\t", 
              signif(loss, 2))
  }
  return(list(lambda = lambda, beta = beta))
}

#####
### Old versions
#####


## Generic Gauss-Newton coordinate descent algorithm


GNcoordinateDescent <- function(beta, y, xi, dxi, 
                                rho = 0.9,
                                c = 1e-4,
                                lambda = 1,
                                reltol = sqrt(.Machine$double.eps),
                                trace = 0) {
  if(trace > 0)
    message("lambda\tdf\tpenalty\t\tloss+penalty")
  beta <- matrix(beta, nrow = length(beta), ncol = length(lambda))
  pIndex <- seq_len(nrow(beta))
  deltak <- 1
  lambda <- sort(lambda, decreasing = TRUE)
  for(j in seq_along(lambda)){
    if(j > 1) 
      beta[, j] <- beta[, j - 1]
    pen <- lambda[j] * sum(abs(beta[, j]))
    r <- xi(beta[, j]) - y
    ell <-  drop(crossprod(r)) 
    loss <- ell + pen
    repeat {
      oldloss <- loss
      for(i in pIndex) {
        d <- dxi(beta[, j], i)
        gamma <- 2 * drop(crossprod(r, d))
        if(beta[i, j] == 0) {
          if(abs(gamma) > lambda[j]) {
            alpha <- drop(crossprod(d, d))
            if(gamma < - lambda[j]) {
              delta <- - (lambda[j] + gamma) / (2 * alpha)
              cg <- c * (gamma + lambda[j])
            } else {  ## gamma > lambda[j]
              delta <- (lambda[j] - gamma) / (2 * alpha)
              cg <- c * (gamma - lambda[j])
            }
            ## Backtracking
            oldval <- ell
            beta[i, j] <- deltak * delta
            r <- xi(beta[, j]) - y
            ell <- drop(crossprod(r))
            val <-  ell + lambda[j] * abs(beta[i, j])
            while(val > oldval + beta[i, j] * cg) {
              deltak <- rho * deltak
              beta[i, j] <- deltak * delta
              r <- xi(beta[, j]) - y
              ell <- drop(crossprod(r))
              val <-  ell + lambda[j] * abs(beta[i, j])
            }
            pen <- pen + lambda[j] * abs(beta[i, j])
          }
        } else { ## beta[i, j] != 0
          beta0i <- beta[i, j]
          cg <- c * (gamma + sign(beta0i) * lambda[j])
          alpha <-  drop(crossprod(d, d))
          di <- 2 * alpha * beta0i - gamma
          if(abs(di) <= lambda[j]) {
            delta <- - beta[i, j]
          } else { ## abs(di) > lambda[j]
            if(di < - lambda[j]) {
              delta <- (lambda[j] - gamma) / (2 * alpha)
            } else { ## di > lambda[j]
              delta <- - (lambda[j] + gamma) / (2 * alpha)
            }
          }
          ## Backtracking
          oldval <- ell + lambda[j] * abs(beta[i, j])
          beta[i, j] <- beta0i + deltak * delta
          r <- xi(beta[, j]) - y
          ell <-  drop(crossprod(r))
          val <-  ell + lambda[j] * abs(beta[i, j])
          while(val > oldval + deltak * delta * cg) {
            deltak <- rho * deltak
            beta[i, j] <- beta0i + deltak * delta
            r <- xi(beta[, j]) - y
            ell <-  drop(crossprod(r))
            val <-  ell + lambda[j] * abs(beta[i, j])
          }
          pen <- pen + lambda[j] * (abs(beta[i, j]) - abs(beta0i))
        }
      }
      loss <- ell + pen
      if(trace == 2)
        message(lambda[j], "\t", sum(abs(beta[,j]) > 0), "\t", pen, "\t", loss)
      if(oldloss - loss < reltol * (abs(loss) + reltol))
        break
    }
    if(trace == 1)
      message(lambda[j], "\t", sum(abs(beta[, j]) > 0), "\t", pen, "\t", loss)
  }
  return(list(lambda = lambda, beta = beta))
}



## An attempt to write a better version using cubic approximations.

GNcoordinateDescentMod <- function(beta, f, gr, quad,
                                lambda = 1,
                                reltol = sqrt(.Machine$double.eps),
                                epsilon = 0,
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
        if(i == 20)
          browser()
        if(beta[i, j] == 0) {
          if(abs(grad[i]) > lambda[j]) {
            beta0 <- beta[ , j]
            alpha <-  2*quad(i, beta[ , j]) + epsilon
            if(grad[i] < -lambda[j]) {
              beta[i, j] <- -(lambda[j] + grad[i])/alpha
            } else {  ## grad[i] > lambda[j]
              beta[i, j] <- (lambda[j] - grad[i])/alpha
            }
            pen <- pen + lambda[j]*abs(beta[i, j])
            grad <- gr(beta[ , j])
          }
        } else { ## beta[i, j] != 0
          beta0i <- beta[i, j]
          alpha <-  2*quad(i, beta[ , j]) + epsilon
          di <- grad[i] - alpha*beta0i
          ## Stick to quadratic approximation
          if(abs(di) > lambda[j]) {
            if (di < -lambda[j]) {
              beta[i, j] <- beta0i - (lambda[j] + grad[i])/alpha
            } else { ## di > lambda[j] 
              beta[i, j] <- beta0i + (lambda[j] - grad[i])/alpha
            }
          } else { ## Turn to cubic approximation
            beta0 <- beta[, j]
            beta0[i] <- 0
            di <- gr(beta0)[i]
            aa <- (alpha*beta0i + di - grad[i])/beta0i^2
            bb <- alpha - 2*aa*beta0i
            ## Is 0 a local minima or is di and beta0i of the same sign?
            if(abs(di) <= lambda[j] || di*beta0i > 0) {
              beta[i, j] <- 0
            } else { ## abs(di) > lambda[j] and di and beta0i of opposite sign
              if (di < -lambda[j]) {
                ##   if(grad[i] >= 0) {
                ##     beta[i, j] <-  - beta0i*(lambda[j] + di)/(grad[i] - lambda[j] - di)
                ##   } else { ## grad[i] < 0
                ##     alpha <-  2*quad(i, beta[ , j])
                beta[i, j] <- (-bb - sqrt(bb^2 - 4*aa*(di + lambda[j])))/(2*aa)
                ## beta[i, j] <- beta0i - (lambda[j] + grad[i])/alpha
              } else { ## di > lambda[j] and di and beta0i of opposite sign
                ## if(grad[i] <= 0) {
                ##  beta[i, j] <- beta0i*(lambda[j] - di)/(grad[i] - di + lambda[j])
                ## } else {
                ##  alpha <-  2*quad(i, beta[ , j])
                beta[i, j] <- (-bb + sqrt(bb^2 - 4*aa*(di - lambda[j])))/(2*aa)
                ## beta[i, j] <- beta0i + (lambda[j] - grad[i])/alpha
              }
            }
          }
          pen <- pen + lambda[j]*(abs(beta[i, j]) - abs(beta0i))
          grad <- gr(beta[ , j])
        }
      }
      val <- f(beta[ , j]) + pen
      if(val <= oldval && oldval-val < reltol*(abs(val) + reltol))
        break
      if(trace == 2)
        cat(lambda[j], "\t", sum(abs(beta[,j]) > 0), "\t", pen, "\t", val, "\n")
    }
    if(trace == 1)
      cat(lambda[j], "\t", sum(abs(beta[,j]) > 0), "\t", pen, "\t", val, "\n")
  }
  return(list(lambda = lambda, beta = beta))
}

## Non-optimal solution based on 'optimize'. This solution solves the
## coodinate wise optimizations in a generic way based on the
## 'optimize' function. It is relatively slow.

coordinateDescentSlow <- function(beta, f, gr, lambda = 1, reltol = sqrt(.Machine$double.eps), epsilon = 0, gamma = 1e-2, trace = 0) {
  if(trace > 0)
    cat("lambda\tdf\tpenalty\t\tloss+penalty\n")
  gamma0 <- gamma
  beta <- matrix(beta, nrow = length(beta), ncol = length(lambda))
  pIndex <- seq_len(dim(beta)[1])
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
          if(abs(grad[i]) > lambda[j] + epsilon) {
            beta0 <- beta[ , j]
            if(grad[i] < -lambda[j]) {
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) + lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(0, gamma0))$minimum
            } else {
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
          beta0[i] <- 0
          if (abs(gr(beta0)[i]) <= lambda[j]) {
            beta[i, j] <- 0
          } else {
            if (beta0i > 0) {
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) + lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(0, beta0i + gamma0))$minimum
            } else { ## beta0i < 0
              tmpf <- function(b) {
                beta0[i] <- b
                f(beta0) - lambda[j]*b
              }
              beta[i, j] <- optimize(tmpf, c(beta0i - gamma0, 0))$minimum
            }
          }
          pen <- pen + lambda[j]*(abs(beta[i, j]) - abs(beta0i))
          val <- f(beta[ , j]) + pen
          grad <- gr(beta[ , j])
        }
      }
      if(val <= oldval && oldval-val < reltol*(abs(val) + reltol))
        break
      if(trace == 1)
        cat(lambda[j], "\t", sum(abs(beta[,j]) > 0), "\t", pen,  "\t", val, "\n")    }
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
