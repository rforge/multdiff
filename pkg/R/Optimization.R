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

## Misc. helper functions. Not exported or documented. ----

penalty <- function(beta, lambda) 
  sum(lambda * abs(beta))

## coordinateDescentQuad ----
##' l1-penalized coordinate descent 
##'
##' Optimizes an l1-penalized loss using a coordinate wise descent algorithm.
##' 
##' The function computes a matrix of optimal parameter values. Each column corresponds to a
##' value of the penalty parameter in \code{lambda}. The estimates are computed in decreasing order 
##' of the penalty parameters, and for each column the previous is used as a warm start.
##' 
##' The algorithm relies on iterative optimization of an l1-penalized 
##' quadratic approximation of the loss using a standard coordinate wise descent algorithm.  
##' A coordinate wise backtracking step is added to ensure that the algorithm takes descent steps. 
##' 
##' This function relies on three auxiliary functions. The loss function \code{f}, its gradient \code{gr}
##' and a third function, \code{quad}, that computes the coefficient of the quadratic approximation 
##' for the coordinate wise optimization. 
##' 
##' The function returns a list with the vector of lambda values as the first entry 
##' and the estimated beta parameters as a matrix in the second entry. Each column 
##' in the matrix corresponds to one lambda value.
##' 
##' @param f a \code{function}. The loss function. A function of \code{beta}.
##' @param gr a \code{function}. The gradient of the loss. A function of \code{beta}.
##' @param quad a \code{function}. The coordinate wise quadratic approximation term. 
##'         A function of the coordinate index and \code{beta}.
##' @param p a \code{numeric}. The length of the parameter vector.        
##' @param beta a \code{numeric}. The vector of initial parameter values. If \code{p} is missing
##'         \code{beta} must be specified. When \code{p} is specified, the 
##'         value of the \code{beta} argument is ignored, and \code{beta = rep(0, p)}. 
##' @param lambda a \code{numeric}. A sequence of penalties. The default value \code{NULL} 
##'         implies an automatic computation of a suitable sequence.
##' @param nlambda a \code{numeric}. The number of lambda values.
##' @param lambda.min.ratio a \code{numeric}. Determines the smallest value of lambda relative to the largest 
##'         (automatically) determined value of lambda. 
##' @param penalty.factor a \code{numeric}. A vector of penalty weight factors.
##' @param rho a \code{numeric}. Step length control in the backtracking.
##' @param c a \code{numeric}. Sufficient decrease control in the backtracking.
##' @param reltol a \code{numeric}. Controls the convergence criterion.
##' @param trace a \code{numeric}. Values above 0 prints increasing amounts of trace information.
##' @return A \code{list} of length 2. The first entry contains \code{lambda} and the second the
##'         matrix of parameter estimates. Each column in the matrix corresponds to an entry in \code{lambda}.
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export

coordinateDescentQuad <- function(f, 
                                  gr, 
                                  quad,
                                  p,
                                  beta,
                                  lambda = NULL,
                                  nlambda = 50,
                                  lambda.min.ratio = 0.001,
                                  penalty.factor = rep(1, p),
                                  rho = 0.9,
                                  c = 1e-4,
                                  reltol = sqrt(.Machine$double.eps),
                                  trace = 0
                                  ) {
  grtol <- reltol ## TODO: Does grtol need to be a different numeric value?
  if(trace > 0)
    message("lambda\tdf\tpenalty\t\tloss+penalty")
  if(missing(p)) {
    if(missing(beta))
      stop("Either the dimension, 'p', or the vector, 'beta', of starting values must be specified.")
    p <- length(beta)
  } else {
    beta <- rep(0, times = p)
  }
  if (is.null(lambda)) {
    lambda <- max(abs(grad(beta)))
    lambda <- lambda * exp(seq(0, log(lambda.min.ratio), length.out = nlambda))
  } else {
    lambda <- sort(lambda, decreasing = TRUE)
  }
  beta <- matrix(beta, nrow = p, ncol = length(lambda))
  pIndex <- seq_len(nrow(beta))
  for(j in seq_along(lambda)){
    if(j > 1) 
      beta[, j] <- beta[, j - 1]
    lamb <- lambda[j] * penalty.factor
    pen <- penalty(beta[, j], lamb)
    val <-  f(beta[, j]) + pen
    grad = gr(beta[, j])
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

## coordinateDescent ----
##' l1-penalized least squares coordinate descent
##'
##' Computes l1-penalized least squares estimates using a coordinate wise descent algorithm of
##' Gauss-Newton quadratic approximations. 
##' 
##' The function computes a matrix of parameter estimates. Each column corresponds to a
##' value of the penalty parameter in \code{lambda}. The estimates are computed in decreasing order 
##' of the penalty parameters, and for each column the previous is used as a warm start.
##' 
##' The algorithm relies on iterative optimization of the l1-penalized Gauss-Newton
##' quadratic approximation using a standard coordinate wise descent algorithm. An 
##' outer backtracking step is added to ensure that the algorithm takes descent steps. 
##' 
##' This function relies on two auxiliary functions. The functions \code{xi} and \code{dxi}, 
##' which are the mean value map and its derivative, respectively. 
##' 
##' The function returns a list with the vector of lambda values as the first entry 
##' and the estimated beta parameters as a matrix in the second entry. Each column 
##' in the matrix corresponds to one lambda value.
##' 
##' @param y a \code{numeric}. The response vector.
##' @param xi a \code{function}. The mean value map. A vector function of \code{beta}.
##' @param dxi a \code{function}. The derivative of the mean value map. A function of \code{beta}
##'         that returns a matrix. 
##' @param p a \code{numeric}. The length of the parameter vector.        
##' @param beta a \code{numeric}. The vector of initial parameter values. If \code{p} is missing
##'         \code{beta} must be specified. When \code{p} is specified, the 
##'         value of the \code{beta} argument is ignored, and \code{beta = rep(0, p)}. 
##' @param lambda a \code{numeric}. A sequence of penalties. The default value \code{NULL} 
##'         implies an automatic computation of a suitable sequence.
##' @param nlambda a \code{numeric}. The number of lambda values.
##' @param lambda.min.ratio a \code{numeric}. Determines the smallest value of lambda relative to the largest 
##'         (automatically) determined value of lambda. 
##' @param penalty.factor a \code{numeric}. A vector of penalty weight factors.
##' @param rho a \code{numeric}. Step length control in the backtracking.
##' @param c a \code{numeric}. Sufficient decrease control in the backtracking.
##' @param reltol a \code{numeric}. Controls the convergence criterion.
##' @param trace a \code{numeric}. Values above 0 prints increasing amounts of trace information.
##' @param N a \code{numeric}. The maximal number of interations in the loops. 
##' @return A \code{list} of length 2. The first entry contains \code{lambda} and the second the
##'         matrix of parameter estimates. Each column in the matrix corresponds to an entry in \code{lambda}.
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export

coordinateDescent <- function(y, 
                              xi, 
                              dxi, 
                              p,
                              beta,
                              lambda = NULL,
                              nlambda = 50,
                              lambda.min.ratio = 0.001,
                              penalty.factor = rep(1, p),
                              rho = 0.9,
                              c = 1e-4,
                              reltol = sqrt(.Machine$double.eps),
                              trace = 0, 
                              N = 1000) {
  kk <- k <- 1
  eps <- .Machine$double.eps
  if(trace > 0)
    message("lambda\tdf\tpenalty\t\tloss+penalty")
  if(missing(p)) {
    if(missing(beta))
      stop("Either the dimension, 'p', or the vector, 'beta', of starting values must be specified.")
    p <- length(beta)
  } else {
    beta <- rep(0, times = p)
  }
  if (is.null(lambda)) {
    grad <- 2 * crossprod(dxi(beta), xi(beta) - y)
    lambda <- max(abs(grad))
    lambda <- lambda * exp(seq(0, log(lambda.min.ratio), length.out = nlambda))
  } else {
    lambda <- sort(lambda, decreasing = TRUE)
  }
  beta <- matrix(beta, nrow = p, ncol = length(lambda))
  pIndex <- seq_len(p)
  delta <- numeric(p)
  beta0 <- numeric(p)
  deltak <- 1
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

## coordinateDescentMF ----
##' l1-penalized matrix free least squares coordinate descent
##'
##' Computes l1-penalized least squares estimates using a coordinate wise descent algorithm of
##' Gauss-Newton quadratic approximations using a matrix free algorithm.
##' 
##' The function computes a matrix of parameter estimates. Each column corresponds to a
##' value of the penalty parameter in \code{lambda}. The estimates are computed in decreasing order 
##' of the penalty parameters, and for each column the previous is used as a warm start.
##' 
##' The algorithm relies on iterative optimization of the l1-penalized Gauss-Newton
##' quadratic approximation using a standard coordinate wise descent algorithm. An 
##' outer backtracking step is added to ensure that the algorithm takes descent steps. This
##' algorithm is matrix free, which means that it does not internally operate with the entire 
##' derivative of the mean value map. Instead, it relies on auxiliary functions to compute the
##' loss, the gradient of the loss and inner products of columns in the derivative of the mean
##' value map. 
##' 
##' This function relies on three auxiliary functions. The loss function \code{f}, its gradient \code{gr}
##' and a third function, \code{cp}, that computes the cross product of the derivative of the mean value map 
##' one column at a time. 
##' 
##' The function returns a list with the vector of lambda values as the first entry 
##' and the estimated beta parameters as a matrix in the second entry. Each column 
##' in the matrix corresponds to one lambda value.
##' 
##' @param f a \code{function}. The loss function. A function of \code{beta}.
##' @param gr a \code{function}. The gradient of the loss. A function of \code{beta}.
##' @param cp a \code{function}. Computes one column in the cross product of the derivative 
##'         of the mean value map. A function of the coordinate index and \code{beta}.
##' @param p a \code{numeric}. The length of the parameter vector.        
##' @param beta a \code{numeric}. The vector of initial parameter values. If \code{p} is missing
##'         \code{beta} must be specified. When \code{p} is specified, the 
##'         value of the \code{beta} argument is ignored, and \code{beta = rep(0, p)}. 
##' @param lambda a \code{numeric}. A sequence of penalties. The default value \code{NULL} 
##'         implies an automatic computation of a suitable sequence.
##' @param nlambda a \code{numeric}. The number of lambda values.
##' @param lambda.min.ratio a \code{numeric}. Determines the smallest value of lambda relative to the largest 
##'         (automatically) determined value of lambda. 
##' @param penalty.factor a \code{numeric}. A vector of penalty weight factors.
##' @param rho a \code{numeric}. Step length control in the backtracking.
##' @param c a \code{numeric}. Sufficient decrease control in the backtracking.
##' @param reltol a \code{numeric}. Controls the convergence criterion.
##' @param trace a \code{numeric}. Values above 0 prints increasing amounts of trace information.
##' @param N a \code{numeric}. The maximal number of interations in the loops. 
##' @return A \code{list} of length 2. The first entry contains \code{lambda} and the second the
##'         matrix of parameter estimates. Each column in the matrix corresponds to an entry in \code{lambda}.
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export

coordinateDescentMF <- function(f, 
                                gr, 
                                cp, 
                                p,
                                beta,
                                lambda = NULL,
                                nlambda = 50,
                                lambda.min.ratio = 0.001,
                                penalty.factor = rep(1, p),
                                rho = 0.9,
                                c = 1e-4,
                                reltol = sqrt(.Machine$double.eps),
                                trace = 0, 
                                N = 1000) {
  kk <- k <- 1
  eps <- .Machine$double.eps
  if(trace > 0)
    message("lambda\tdf\tpenalty\t\tloss+penalty")
  if(missing(p)) {
    if(missing(beta))
      stop("Either the dimension, 'p', or the vector, 'beta', of starting values must be specified.")
    p <- length(beta)
  } else {
    beta <- rep(0, times = p)
  }
  if (is.null(lambda)) {
    lambda <- max(abs(gr(beta)))
    lambda <- lambda * exp(seq(0, log(lambda.min.ratio), length.out = nlambda))
  } else {
    lambda <- sort(lambda, decreasing = TRUE)
  }
  beta <- matrix(beta, nrow = p, ncol = length(lambda))
  pIndex <- seq_len(p)
  delta <- numeric(p)
  beta0 <- numeric(p)
  deltak <- 1
  for(j in seq_along(lambda)){
    if(j > 1) 
      beta[, j] <- beta[, j - 1]
    lamb <- lambda[j] * penalty.factor
    pen <- penalty(beta[, j], lamb)
    ell <-  f(beta[, j]) 
    loss <- ell + pen
    repeat {
      oldloss <- loss
      beta0 <- beta[, j]
      delta <- numeric(p)
      cpX <- vector("list", p)  ## List of columns of cross products
      dnonzero <- logical(p)
      gamma0 <- gamma <- gr(beta0) 
      sqloss <- penalty(beta0, lamb) 
      repeat {
        sqoldloss <- sqloss
        for (i in pIndex) {
          if (beta0[i] == 0) {
            if (abs(gamma[i]) > lamb[i]) {
              if (is.null(cpX[[i]])) {
                dnonzero[i] <- TRUE
                cpX[[i]] <- as.vector(cp(i, beta0))
              }
              alpha <- cpX[[i]][i] 
              if (gamma[i] < - lamb[i]) {
                d <- - (lamb[i] + gamma[i]) / (2 * alpha)
              } else {  ## gamma > lambda[j]
                d <- (lamb[i] - gamma[i]) / (2 * alpha)
              }
              gamma <- gamma + (2 * d) * cpX[[i]]
              delta[i] <- delta[i] + d
              beta0[i] <- beta0[i] + d
            }
          } else { ## beta0[i] != 0 
            if (is.null(cpX[[i]])) {
              dnonzero[i] <- TRUE
              cpX[[i]] <- as.vector(cp(i, beta0))
            }
            alpha <- cpX[[i]][i] 
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
            gamma <- gamma + (2 * d) * cpX[[i]]
            delta[i] <- delta[i] + d
            beta0[i] <- beta0[i] + d
          }
        }
        pen <- penalty(beta0, lamb) 
        deltanz <- delta[dnonzero]
        cpXd <- sapply(cpX[dnonzero], function(a) crossprod(a[dnonzero], deltanz))
        if (length(cpXd) > 0) {
        sqloss <- drop(crossprod(deltanz, cpXd) + 
                         crossprod(gamma[dnonzero], deltanz)) + pen
        } 
        if (trace >= 4)
          message("Quadratic loop:", 
                  signif(lambda[j], 2), 
                  "\t", sum(abs(beta0) > 0), 
                  "\t", signif(pen, 2), 
                  "\t", signif(sqloss, 2))
        if (sqoldloss - sqloss < reltol * (abs(sqloss) + reltol))
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
      cg <- c * (crossprod(gamma0, delta) + pen - oldpen)
      beta0 <- beta[, j] + deltak * delta
      ell <- f(beta0)
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
        ell <-  f(beta0)
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
