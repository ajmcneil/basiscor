# FUNCTIONS FOR POPULATION BASIS CORRELATIONS

#' Compute Basis Correlation
#'
#' @param object an object of class copula or Copula, or a data matrix, or a bivariate function describing a copula.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param ... other arguments to function.
#'
#' @return value of polynomial rank correlation.
#' @export
#'
#' @examples
#' basiscor(copula::claytonCopula(2), 2, 2)
basiscor <- function(object, j = 1L, k = 1L, ...){
  if (methods::is(object, "matrix"))
    basiscordata(object, j, k, ...)
  else if ((methods::is(object, "copula")) | (methods::is(object, "Copula")))
    basiscorcopula(object, j, k, copobj = TRUE, ...)
  else if (methods::is(object, "function"))
    basiscorcopula(object, j, k, copobj = FALSE, ...)
  else
    stop("Unknown object type")
}

#' Compute Basis Correlation for Copula Object or Function
#'
#' @param copula an object of class parCopula or a function.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param copobj logical parameter for copula object.
#' @param method method of calculation which can be "p" or "d".
#' @param type type of basis correlation can be legendre or cosine.
#' @param ... other arguments to function.
#'
#' @return value of polynomial rank correlation.
#' @export
#' @import copula
#'
#' @examples
#' basiscorcopula(copula::claytonCopula(2), 2, 2, copobj = TRUE)
basiscorcopula <- function(copula, j = 1L, k = 1L, copobj, method = "p", 
                          type = "legendre", ...){
  j <- as.integer(j)
  k <- as.integer(k)
  if ((j == 0) & (k == 0))
    return(1)
  if ((j == 0) | (k == 0))
    return(0)
  if ((j < 0) | (k < 0) | (j%%1 != 0) | (k%%1 != 0))
    stop("j and k must be positive integers")
  if (copobj){
    if (!methods::is(copula, "parCopula"))
    stop("Function requires a parametric copula object")
    if (methods::is(copula, "tCopula"))
    if (copula@parameters[2]%%1 != 0) #pCopula not implemented for non-integer df
      method <- "d"
  }
  result <- stats::integrate(basiscor_outer, lower = 0, upper = 1, j = j, k = k, copobj = copobj,
                      copula = copula, method = method, type = type, ...)$value
  if (method == "p"){
    if (type == "cosine")
      result <- result -2
    if (type == "legendre"){
      scale <- j * (j+1) * k * (k+1)
      result <-  (result * scale - 1) * sqrt(2 * j + 1) * sqrt(2 * k + 1)
    }
  }
  result
}

#' Outer Integrand for Basis Correlation
#'
#' @param v vector argument of function.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param copobj logical parameter for copula object.
#' @param copula an object of class parCopula or a function.
#' @param method method of calculation which can be "p" or "d".
#' @param type type of basis correlation can be legendre or cosine.
#'
#' @return value of outer integrand
#' @keywords internal
#'
basiscor_outer <- function(v, j, k, copobj, copula, method, type, ...){
  out <- rep(NA, length(v))
  for (i in 1:length(v))
  {
    tmp <- stats::integrate(basiscor_inner, lower = 0, upper = 1, 
                     v = v[i], j = j, k = k, copobj = copobj, copula = copula, 
                     method = method, type = type, ...)
    out[i] <- tmp$value
  }
  out
}

#' Inner Integrand for Basis Correlation
#'
#' @param u vector argument of function.
#' @param v vector argument of function.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param copobj logical parameter for copula object.
#' @param copula an object of class parCopula or a function.
#' @param method method of calculation which can be "p" or "d".
#' @param type type of basis correlation can be legendre or cosine.
#'
#' @return value of inner integrand
#' @keywords internal
basiscor_inner <- function(u, v, j, k, copobj, copula, method, type, ...){
  if ((copobj) & (method == "p"))
    part1 <- pCopula(cbind(u, v), copula)
  else if ((copobj) & (method == "d"))
    part1 <- dCopula(cbind(u, v), copula)
  else
    part1 <- copula(u, v, ...) 
  if (method == "d")
    output <- switch(type,
                     cosine = 2 * cos(j * pi * u) * cos(k * pi * v) * (-1)^(j + k) *  part1,
                     sLegendre(u, j) * sLegendre(v, k) * sqrt(2 * j + 1) * sqrt(2 * k + 1) * part1
    )
  if (method == "p")
    output <- switch(type,
         "cosine" = 2 * j * k * (pi^2) * sin(j * pi * u) * sin(k * pi * v) * (-1)^(j + k) * part1,
        "legendre" = sLegendre(u,j, deriv = TRUE) * sLegendre(v, k, deriv = TRUE) * part1 / (j * (j + 1) * k * (k + 1))
        )
  output
}

#' Compute Matrix of Basis Correlations
#'
#' @param object an object of class copula or Copula, or a data matrix, or a bivariate function describing a copula.
#' @param maxorder maximum order of the polynomials.
#' @param symmetric logical variable stating whether matrix is known a priori to be symmetric.
#' @param ... other parameters passed to underlying functions.
#'
#' @return a square matrix with number of rows and columns equal to maxorder
#' @export
#'
#' @examples
#' basiscormatrix(copula::claytonCopula(2))
basiscormatrix <- function(object, maxorder = 4, symmetric = FALSE, ...){
  output <- matrix(NA, nrow = maxorder, ncol = maxorder)
  for (j in 1:maxorder){
    k1 <- 1
    if (symmetric) 
      k1 <- j
    for (k in k1:maxorder)
      output[j,k] <- basiscor(object, j, k, ...)
  }
  if (symmetric)
    output[lower.tri(output)] <- t(output)[lower.tri(t(output))]
  output
}

#' Extremal Basis Correlation for Two Shifted Legendre Polynomials
#'
#' The maximum (`case = "max"`) or minimum (`case = "min"`) value attainable
#' by the population basis correlation \code{basiscor(copula, j, k)} over all
#' copulas, given by matching the quantile functions of \eqn{P_j(U)} and
#' \eqn{P_k(U)} comonotonically (max) or countermonotonically (min), where
#' \eqn{P_j = \sqrt{2j+1} L_j} is the orthonormal shifted Legendre polynomial
#' of degree \code{j}. Uses \code{udp::udpquantile()}, the quantile function
#' of \eqn{L_j(U)} carried by \code{udp::udplegendre(j)}.
#'
#' @param j degree of first polynomial
#' @param k degree of second polynomial
#' @param case character variable which should be "max" or "min"
#'
#' @return value of extremal Legendre correlation
#' @export
#'
#' @examples
#' extremalLegendre(3, 4)
extremalLegendre <- function(j, k, case = "max") {
  Tj <- udp::udplegendre(j)
  Tk <- udp::udplegendre(k)
  integrand <- switch(case,
    max = function(u) udp::udpquantile(Tj, u) * udp::udpquantile(Tk, u),
    min = function(u) udp::udpquantile(Tj, u) * udp::udpquantile(Tk, 1 - u),
    stop("Unknown case for extremum")
  )
  sqrt(2 * j + 1) * sqrt(2 * k + 1) * stats::integrate(integrand, 0, 1)$value
}

#' Compute Sample Basis Correlation
#'
#' @param data a matrix of data wit two columns.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#' @param method method of calculation
#'
#' @return sample polynomial rank correlation value.
#' @export
#'
basiscordata <- function(data, j, k, type = "legendre", method = "T1"){
  if (! methods::is(data, "matrix"))
    stop("Must supply data matrix")
  if (ncol(data) != 2)
    stop("Data matrix should have 2 columns exactly")
  R1 <- rank(data[,1])
  R2 <- rank(data[,2])
  if ((type != "legendre") & (method == "T6"))
    stop("Discrete polynomials only available for Legendre correlations")
  n <- length(R1)
  switch(method,
         T1 = mean(basisfunc(R1/(n+1), j, type) * basisfunc(R2/(n+1), k, type)),
         T2 = mean(basisfunc((R1-0.5)/n, j, type) * basisfunc((R2-0.5)/n, k, type)),
         T3 = stats::cor(basisfunc(R1/(n+1), j, type), basisfunc(R2/(n+1), k, type)),
         T4 = stats::cor(basisfunc((R1-0.5)/n, j, type), basisfunc((R2-0.5)/n, k, type)),
         T5 = n * sum((basisintegral(R1/n, j, type) - basisintegral((R1-1)/n, j, type))*
                      (basisintegral(R2/n, k, type) - basisintegral((R2-1)/n, k, type))),
         T6 = sqrt((2*j+1)*(2*k+1))*exp((lff(n-1,j)+lff(n-1,k)-lff(n+j,j)-lff(n+k,k))/2) *
           mean(discreteLegendre(R1,n,j)*discreteLegendre(R2,n,k)))
}

#' Log of Falling Factorial Function
#'
#' @param r first parameter
#' @param k second parameter
#'
#' @return value of function
#' @keywords internal
#'
lff <- function(r,k){
  lgamma(r+1)-lgamma(r-k+1)
}

#' Discrete Legendre Polynomial
#' 
#'
#' @param r rank of observation in sample.
#' @param n size of sample/
#' @param degree non-negative integer giving degree of polynomial.
#'
#' @return vector of values of polynomial.
#' @export
#'
#' @examples
#' discreteLegendre(7, 20, 2)
discreteLegendre <- function(r, n, degree){
  output <- 0
  for (k in 0:degree){
    output <- output + choose(degree,k)* choose(degree+k,k) * exp(lff(r-1,k) -lff(n-1,k)) * (-1)^k
  }
  output*(-1)^degree
}