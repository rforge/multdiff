### R code for construction and computation with multivariate process data.
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

## Process constructor ----
##' Multivariate process models constructor
##' 
##' Constructs multivariate models for dynamic processes.
##' 
##' @param y a \code{numeric}. A response matrix of dimensions d by m.
##' @param x a \code{numeric}. A predictor matrix of dimensions d by m.
##' @param t a \code{numeric}. A single time step size, or a vector of time step sizes of length m. 
##' @param suff a \code{logical}. Should sufficient statistics and the MLE be computed.
##' @param weights a \code{character} or \code{numeric}. Specifies parameter weights. Default value 
##'        \code{'unit'} means that all parameters are weighted equally with a unit weight, 
##'        whereas \code{'adaptive'} means that the parameters are weighted inverse proportionally 
##'        to the MLE (will only be used if \code{suff = TRUE}). Alternatively, a numeric specifying 
##'        the individual parameter weights.
##' @return an object of class \code{multModel} containing the elements
##'  \tabular{ll}{
##'  y \tab the d by m data matrix \cr
##'  x \tab the d by m data matrix \cr
##'  t \tab the time step(s) \cr
##'  ss \tab sum of squares (computed if \code{suff = TRUE}) \cr
##'  xyt \tab d by d crossproduct of x and y (computed if \code{suff = TRUE}) \cr
##'  xxt \tab d by d crossproduct of x (computed if \code{suff = TRUE}) \cr
##'  B \tab the MLE if \code{suff = TRUE}, and the zero matrix otherwise \cr
##'  w \tab a matrix of parameter weights \cr
##'  d \tab dimension of the process \cr
##'  m \tab number of d-dimensional observations \cr
##'  n \tab equals \code{m * d} \cr
##'  p \tab quals \code{d * d} \cr
##'  sigma \tab the observation variance parameter \cr
##'  status \tab 0 means no problems or errors encountered \cr
##'  msg \tab error messages or other notable conditions 
##'  }
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export

multModel <- function(y, x, t, suff = TRUE, weights = "unit", sigma = 1) {
  if (!(is.matrix(y) & is.matrix(x) & all(dim(x) == dim(y))))
    stop("Arguments 'y' and 'x' must be matrices of the same dimensions.")
  m <- ncol(y)
  d <- nrow(y)
  n <- m * d
  p <- d * d
  
  if (length(t) > 1)
    stop("Multiple time steps is currently not supported.")
  
  if (is.numeric(weights)) {
    if (length(weights) != d * d)
      warning("The length of 'weights' does not match the number of parameters. The weights are recycled.")
    w <- matrix(weights, d, d)
  } else {
    w <- matrix(1, d, d)
  }
  
  B <- matrix(0, d, d) 
  status <- 0
  msg <- ""
  
  if (suff) {
    ss <- sum(y^2)
    xyt <- tcrossprod(x, y)
    xxt <- tcrossprod(x)
    MLEtry <- try(solve(xxt, xyt), silent = TRUE)
    if (inherits(MLEtry, "try-error")) {
      status <- 2
      msg <- MLEtry
    } else {
      Btry <- try(logm(t(MLEtry)) / t, silent = TRUE)
      if (inherits(Btry, "try-error")) {
        status <- 2
        msg <- Btry
      } else {
        B <- Btry
        if (weights[1] == 'adaptive')
          w <- 1 / (abs(B) + 1e-8)
      }
    }
  } else {  ## Don't compute sufficient stats
    ss <- xyt <- xxt <- 0 
  }
    
    structure(list(y = y,
                   x = x,
                   t = t,
                   ss = ss,
                   xyt = xyt,
                   xxt = xxt,
                   B = B,
                   w = w,
                   d = d,
                   m = m,
                   n = n,
                   p = p,
                   suff = suff,
                   dfMethod = "pen",
                   sigma = sigma,
                   status = status,
                   msg = msg),
              class = "multModel")
}

## Misc. functions that are not exported ----

xi <- function(model, B) {
  d <- model$d
  t <- model$t
  x <- model$x
  dim(B) <- c(d, d)
  as.numeric(expm(t * B) %*% x)
}

dxi <- function(model, B) {
  d <- model$d
  t <- model$t
  x <- model$x
  p <- model$p
  n <- model$n
  dim(B) <- c(d, d)
  F <- numeric(p)
  dim(F) <- dim(B)
  dxi <- matrix(0, n, p)
  for(j in 1:p) {
    F[j] <- 1
    dxi[, j] <- as.vector(expmFrechet(t * B, F, expm = FALSE)$Lexpm %*% x)
    F[j] <- 0
  }
  t * dxi
}

crossprodCol <- function(model, i, B) {
  d <- model$d
  t <- model$t
  dim(B) <- c(d, d)
  if (model$suff) {
    xxt <- model$xxt
    E <- matrix(0, d, d)
    E[i] <- 1
    LB <- expmFrechet(t * B, E, expm = FALSE)$Lexpm
    as.vector(t^2 * expmFrechet(t * t(B), LB %*% xxt, expm = FALSE)$Lexpm)
  } else {
    X <- dxi(model, B) 
    drop(crossprod(X, X[, i]))
  }
}


quad <- function(model, i, B) {
  d <- model$d
  t <- model$t
  dim(B) <- c(d, d)
  if (model$suff) { 
    xxt <- model$xxt
    E <- matrix(0, d, d)
    E[i] <- 1
    LB <- expmFrechet(t * B, E, expm = FALSE)$Lexpm
    t^2 * sum(LB * (LB %*% xxt))
  } else {
    p <- model$p
    x <- model$x
    F <- numeric(p)
    dim(F) <- dim(B)
    F[i] <- 1
    dxi <- as.vector(expmFrechet(t * B, F, expm = FALSE)$Lexpm %*% x)
    t^2 * sum(dxi^2)
  }
}

dfPen <- function(model, B, ...) {
  d <- model$d
  t <- model$t
  y <- model$y
  x <- model$x
  n <- model$n
  dim(B) <- c(d, d)
  tmp <- matrix(0, 3 * d, 3 * d)
  tmp[1:d, 1:d] <- tmp[1:d + d, 1:d + d] <- 
    tmp[1:d + 2 * d, 1:d + 2 * d] <- t * B
  r <- y - expm(t * B) %*% x
  xr <- x %*% t(r)
  ii <- which(B != 0)
  if (length(ii) > 0) {
    pp <- length(ii)
    Dzeta <- matrix(0, n, pp)
    J <- matrix(0, pp, pp)
    E1 <- matrix(0, d, d)
    for(i in seq_along(ii)) {
      E1[ii[i]] <- 1
      Dzeta[, i] <- t * expmFrechet(t * B, E1, expm = FALSE)$Lexpm %*% x
      tmp[1:d, 1:d + d] <- E1
      tmp[1:d + d, 1:d + 2 * d] <- xr
      J[i, ] <- - t^2 * as.vector(t(expm(tmp)[1:d, 1:d + 2 * d]))[ii]
      E1[ii[i]] <- 0
    }
    J <- J + t(J)
    G <- crossprod(Dzeta)
    J <- G + J
    JinvG <- try(solve(J, G), silent = TRUE)
    if (inherits(JinvG, "try-error")) {
      warning("Error in matrix inversion when computing the effective degrees of freedom. Number of non-zero parameters returned.")
      df <- length(ii)
    } else {
      df <- sum(diag(JinvG))
    }
  } else {
    df <- 0
  }
  df
}

dfCon <- function(model, B, omega = rep(1, model$p), ...) {
  d <- model$d
  t <- model$t
  y <- model$y
  x <- model$x
  n <- model$n
  dim(B) <- c(d, d)
  tmp <- matrix(0, 3 * d, 3 * d)
  tmp[1:d, 1:d] <- tmp[1:d + d, 1:d + d] <- 
    tmp[1:d + 2 * d, 1:d + 2 * d] <- t * B
  r <- y - expm(t * B) %*% x
  xr <- x %*% t(r)
  ii <- which(B != 0)
  if (length(ii) > 0) {
    pp <- length(ii)
    Dzeta <- matrix(0, n, pp)
    J <- matrix(0, pp, pp)
    E1 <- matrix(0, d, d)
    for(i in seq_along(ii)) {
      E1[ii[i]] <- 1
      Dzeta[, i] <- t * expmFrechet(t * B, E1, expm = FALSE)$Lexpm %*% x
      tmp[1:d, 1:d + d] <- E1
      tmp[1:d + d, 1:d + 2 * d] <- xr
      J[i, ] <- - t^2 * as.vector(t(expm(tmp)[1:d, 1:d + 2 * d]))[ii]
      E1[ii[i]] <- 0
    }
    J <- J + t(J)
    G <- crossprod(Dzeta)
    J <- G + J
    gamma <- as.vector(omega[ii] * sign(B[ii]))
    JinvG <- try(solve(J, G), silent = TRUE)
    JinvGam <- try(solve(J, gamma), silent = TRUE)
    if (inherits(JinvG, "try-error") || inherits(JinvGam, "try-error")) {
      warning("Error in matrix inversion when computing the effective degrees of freedom. Number of non-zero parameters minus 1 returned.")
      df <- length(ii) - 1
    } else {
      df <- sum(diag(JinvG)) - sum(gamma * (JinvG %*% JinvGam)) / (sum(gamma * JinvGam))
    }
  } else {
    df <- 0
  }
  df
}

## Misc. generic functions and methods ----

##' @export
dim.multModel <- function(x)
  dim(x$y)

##' Returns the response 
##' 
##' @param object an object from which the response is to be obtained.
##' @seealso \code{\link{multModel}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
response <- function(model, ...)
  UseMethod("response")

##' @rdname response
##' @export
response.multModel <- function(model, ...)
  model$y

##' Prediction for multivariate models
##' 
##' @param object an object of class \code{multModel}.
##' @param par the parameter vector.
##' @seealso \code{\link{predict}}, \code{\link{multModel}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @importFrom stats predict
##' @export
predict.multModel <- function(object, par, ...) 
  xi(object, par)


##' Derivative of the predictor 
##' 
##' @param object an object of class \code{multModel}.
##' @param par the parameter vector.
##' @return An n by p matrix.
##' @seealso \code{\link{predict.multModel}}, \code{\link{multModel}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
dpredict <- function(object, par, ...)
  UseMethod("dpredict")

##' @rdname dpredict
##' @export
dpredict.multModel <- function(object, par, ...) 
  dxi(object, par)

##' Computation of the squared error loss function
##' 
##' @param object for which the loss is to be computed.
##' @param par the parameter vector.
##' @return a numeric.
##' @seealso \code{\link{grad}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
loss <- function(model, par, ...)
  UseMethod("loss")

##' @export
loss.default <- function(model, par, ...) {
  y <- response(model)
  pred <- predict(model, par)
  sum((y - pred)^2)
}

##' @rdname loss
##' @export
loss.multModel <- function(model, par = model$B, ...) {
  d <- model$d
  t <- model$t
  B <- par
  dim(B) <- c(d, d)
  if (model$suff) {
    ss <- model$ss
    xyt <- model$xyt
    xxt <- model$xxt
    eB <- expm(t * B)
    loss <- ss - 2 * sum(t(eB) * xyt) + sum(eB * (eB %*% xxt))
  } else {
    y <- model$y
    x <- model$x
    loss <- sum((y - expm(t * B) %*% x)^2)
  }
  loss
}

##' Computation of the gradient of the loss function
##' 
##' @param model for which the gradient of the loss is to be computed.
##' @param par the parameter vector.
##' @param ...
##' @return a numeric vector.
##' @seealso \code{\link{loss}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
grad <- function(model, par, ...)
  UseMethod("grad")

##' @export
grad.default <- function(model, par, ...) {
  r <- as.vector(predict(model, par) - response(model))
  2 * drop(r %*% dpredict(model, par))
}

##' @rdname grad
##' @export
grad.multModel <- function(model, par = model$B, ...) {
  d <- model$d
  t <- model$t
  dim(par) <- c(d, d)
  if (model$suff) {
    ss <- model$ss
    xyt <- model$xyt
    xxt <- model$xxt
    F <- xyt - xxt %*% t(expm(t * par))
    grad <- as.vector(- 2 * t * t(expmFrechet(t * par, F, expm = FALSE)$Lexpm))
  } else {
    y <- model$y
    r <- as.vector(xi(model, par) - y)
    grad <- 2 * drop(r %*% dxi(model, par))
  }
  grad
}

##' Effective degrees of freedom.
##' 
##' Estimation of the effective degrees of freedom.
##' 
##' @param model the model for which the degrees of freedom is to be computed.
##' @param par a \code{numeric}. The parameter vector.
##' @param method a \code{character}. The method for estimation of degrees of freedom. The default is to use 
##'        the method specified by the model. The possible values are \code{'pen'}, \code{'sub'}, 
##'        \code{'con'} and \code{'nonzero'}.
##' @param ...
##' @return a numeric.
##' @seealso \code{\link{risk}}, \code{\link{riskHat}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
dfEff <- function(model, par, ...)
  UseMethod("dfEff")

##' @export
dfEff.default <- function(model, par, ...)
  length(par)

##' @rdname dfEff
##' @export
dfEff.multModel <- function(model, par = model$B, method = model$dfMethod, ...)
  switch(method,
         pen = dfPen(model, par),
         sub = dfPen(model, par),
         con = dfCon(model, par, model$w),
         nonzero = sum(par != 0),
         NextMethod(model = model, par = par))

##' Computation of the squared error risk.
##' 
##' @param object for which the risk is to be computed.
##' @param par0 the best / true parameter vector.
##' @param par the estimated parameter vector.
##' @param ...
##' @return a numeric.
##' @seealso \code{\link{dfEff}}, \code{\link{riskHat}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
risk <- function(model, par0, par, ...)
  UseMethod("risk")

##' @export
risk.default <- function(model, par0, par, ...) {
  pred0 <- predict(model, par0)  
  pred <- predict(model, par)
  sum((pred - pred0)^2)
}

##' Estimation of the risk.
##' 
##' Computation of an estimate of the quadratic risk based on the effective degrees of 
##' freedom. The estimate can be used for model selection.
##' 
##' @param model the model object for which the risk is to be computed.
##' @param par a \code{numeric}. The estimated parameter vector.
##' @param df a \code{numeric}. The effective degrees of freedom.
##' @param sigma a \code{numeric}. The variance parameter. 
##' @return a numeric.
##' @seealso \code{\link{df}}, \code{\link{risk}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
riskHat <- function(model, par, df, sigma, ...)
  UseMethod("riskHat")

##' @export
riskHat.default <- function(model, par, df = length(par), sigma = 1, ...) 
  loss(model, par) + sigma^2 * (2 * df - model$n)

##' @rdname riskHat
##' @export
riskHat.multModel <- function(model, par = model$B, df = dfEff(model, par), sigma = model$sigma, ...) {
  force(par)
  force(df)
  force(sigma)
  NextMethod(par = par, df = df, sigma = sigma)
  #  t <- object$t
  #  x <- object$x
  #  y <- object$y
  #  d <- nrow(x)
  #  dim(B) <- c(d, d)
  #  sum((y - expm(t * B) %*% x)^2) + sigma^2 * (2 * df - n)
}

##' Fitting a multivariate process model
##' 
##' Estimation of parameters in a multivariate process model. The function fits the models by minimizing 
##' an l1-penalized squared error loss by a coordinate wise descent algorithm. Alternatively, the function 
##' can do unpenalized least squares forward or backward model search.
##' 
##' The stepwise model search algorithms as well as the algorithm for l1-panalized estimation return
##' a sequence of models parametrized by a vector of tuning parameters. For the model search, the tuning 
##' parameter is the dimension of the model (the number of nonzero entries in the parameter vector). For 
##' the l1-panalized models, the tuning parameter is a sequence of penalty parameters. 
##' 
##' The scope argument specifies the smallest and largest model fitted. They must be nested. The zero weights
##' of the \code{multModel} object are used to specify the smallest model, if the scope argument does not 
##' contain a smallest model. 
##' 
##' @param model 
##' @seealso \code{\link{multModel}}, \code{\link{coordinateDescent}}, \code{\link{coordinateDescentMF}},
##'         \code{\link{stepOptim}}.
##' @return a \code{multModel} object with the additional elements
##'  \tabular{ll}{
##'  lambda \tab the sequence of tuning parameters \cr
##'  Blambda \tab the corresponding matrix of estimated parameters. The i'th column corresponds to the i'th tuning parameter. \cr
##'  status \tab either 0 (meaning no errors), 1 (convergence), 2 (MLE not computed) or 3 (other errors) 
##'               indicating different problems or errors \cr
##'  msg \tab errors encountered or other notable conditions
##'  }         
##' @export
fit <- function(model, ...)
  UseMethod("fit")

##' @rdname fit
##' @param method a \code{character}. The default, \code{"cd"}, means a standard coordinate descent algorithm. 
##'         while \code{'mf'} means the same coordinate descent algorithm, but using a matrix free 
##'         implementation. A choice of \code{'forward'} or \code{'backward'} means least squares model search.   
##' @param scope a \code{numeric} or a \code{list}. See \code{\link{stepOptim}}.
##' @param ... additional arguments passed on to the optimization algorithm.  
##' @export        
fit.multModel <- function(model, method = 'cd', scope = seq(1, model$p), ...) {
  if (method %in% c('cd', 'mf', 'forward', 'backward')) {
    if (method == 'cd') 
      result <- try(coordinateDescent(y = as.vector(model$y), 
                                      xi = function(B) xi(model, B),
                                      dxi = function(B) dxi(model, B),
                                      p = model$p,
                                      penalty.factor = as.vector(model$w),
                                      ...), silent = TRUE)
    if (method == 'mf')
      result <- try(coordinateDescentMF(f = function(B) loss(model, B), 
                                        gr = function(B) grad(model, B),
                                        cp = function(i, B) crossprodCol(model, i, B),
                                        p = model$p,
                                        penalty.factor = as.vector(model$w),
                                        ...), silent = TRUE)
    if (method %in% c('forward', 'backward')) {
      B <- numeric(model$p)
      f <- function(par, pattern, ...) {
        B[pattern] <- par
        loss(model, par = B)
      }
      gr <- function(par, pattern, ...) {
        B[pattern] <- par
        grad(model, B)[pattern]  
      }
      if (!is.list(scope)) {
        upper <- scope
        lower <- which(model$w == 0)
        scope <- list(lower = lower, upper = upper)
      }
      result <- try(stepOptim(f = f, 
                              gr = gr, 
                              par = model$B, 
                              direction = method,
                              scope = scope,
                              ...), silent = TRUE)
    }
      
    if (inherits(result, "try-error")) {
      model$status <- 3
      model$msg <- paste(model$msg, result)
    } else {
      model$lambda <- result$lambda
      model$Blambda <- result$beta
      model$status <- max(result$status, model$status)
      if (result$status == 1)
        model$msg <- paste(model$msg, "Coordinate descent algorithm reached max interations for some lambda.")
      if (result$status == 2)
        model$msg <- paste(model$msg, result$msg)
    }
  } else {     
    warning("Method not available. Model returned without fitting.")
  }
  model
}
