# FUNCTIONS FOR POPULATION BASIS CORRELATIONS

# Coerces j, k to integers, validating that both are non-negative whole
# numbers (checked *before* coercion, so a non-integer such as 2.5 is
# rejected rather than silently truncated), and returns the value implied by
# orthonormality when either is 0: 1 if both are 0 (P_0 = 1 is trivially
# self-correlated), 0 if only one is (E[P_0(U) P_k(U)] = 0 for k > 0,
# exactly, for every copula). `shortcut` is NA when the caller must compute
# the general case. Shared by basiscorcopula() and basiscordata() so the two
# agree on this convention.
validate_degrees <- function(j, k) {
  if (!is.numeric(j) || !is.numeric(k) || length(j) != 1 || length(k) != 1 ||
      anyNA(c(j, k)) || j < 0 || k < 0 || j != round(j) || k != round(k)) {
    stop("j and k must be non-negative integers", call. = FALSE)
  }
  j <- as.integer(round(j))
  k <- as.integer(round(k))
  shortcut <- if (j == 0 && k == 0) 1 else if (j == 0 || k == 0) 0 else NA
  list(j = j, k = k, shortcut = shortcut)
}

# n-point Gauss-Legendre quadrature nodes and weights on [0, 1], via the
# Golub-Welsch eigendecomposition of the symmetric tridiagonal Jacobi matrix
# for the Legendre three-term recurrence (off-diagonal beta_k = k / sqrt(4k^2
# - 1); nodes are its eigenvalues, weights the squared first components of
# the corresponding eigenvectors, rescaled from [-1, 1] to [0, 1]).
gauss_legendre_01 <- function(n) {
  k <- seq_len(n - 1)
  beta <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1), 2:n)] <- beta
  J[cbind(2:n, seq_len(n - 1))] <- beta
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  list(nodes = (eig$values[ord] + 1) / 2, weights = eig$vectors[1, ord]^2)
}

# Tensor-product Gauss-Legendre quadrature for basiscorcopula()'s double
# integral, replacing basiscor_outer()'s nested stats::integrate() for both
# method = "p" and method = "d".
#
# method = "p": the integrand is smooth (built from pCopula and the
# derivative of a fixed-degree polynomial/cosine basis function, never the
# copula density). Validated against a 300-node-per-axis reference across
# several copula families at degree up to (7, 8): error <= 6e-5 at ngrid =
# 40, comfortably inside stats::integrate()'s own ~1.2e-4 default tolerance,
# at roughly 1.5-2x the speed -- replacing many small nested adaptive
# integrate() calls with one vectorized evaluation of the copula.
#
# method = "d" uses the copula density directly, which is far more expensive
# per evaluation for elliptical copulas (normal, t) -- exactly the families
# method = "d" exists for, since their pCopula has no closed form -- so the
# same one-vectorized-batch trick is a much bigger win there (~13x at ngrid =
# 40 for a t copula) than for method = "p". It also fixes a real reliability
# problem: stats::integrate() on the nested density integrand returns NA
# outright for some entirely ordinary cases (e.g. normalCopula(0.6) at j = 6,
# k = 6), not just for extreme tail-dependence parameters. But the density
# integrand converges more slowly than "p"'s: for strongly dependent,
# low-degrees-of-freedom t copulas at higher polynomial degree, ngrid = 40
# can be off by ~1% and ngrid = 150-200 is needed to match "p"-level (~1e-4)
# accuracy. basiscorcopula() therefore defaults ngrid to 80 for method = "d"
# (not 40): a compromise that stays reliable and meaningfully faster than the
# old code even in the hard cases, while getting close to "p"-level accuracy
# for the more common, more weakly dependent cases -- see
# basiscorcopula()'s ngrid documentation for how to tighten it further.
basiscorcopula_gl <- function(cop, j, k, copobj, method, type, ngrid, ...) {
  GL <- gauss_legendre_01(ngrid)
  uv <- expand.grid(u = GL$nodes, v = GL$nodes)
  w <- as.vector(outer(GL$weights, GL$weights))
  vals <- basiscor_inner(uv$u, uv$v, j, k, copobj, cop, method, type, ...)
  sum(w * vals)
}

#' Compute Basis Correlation
#'
#' @param object an object of class copula or Copula, or a data matrix, or a bivariate function describing a copula.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param ... further arguments passed on to \code{\link{basiscordata}} (when
#'   \code{object} is a data matrix) or \code{\link{basiscorcopula}}
#'   (otherwise), such as \code{type} (\code{"legendre"} or \code{"cosine"})
#'   and, for \code{basiscordata}, \code{method}, or for
#'   \code{basiscorcopula}, \code{method} (\code{"p"} or \code{"d"}).
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
#' @param cop an object of class parCopula or a function.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param copobj logical parameter for copula object.
#' @param method method of calculation which can be "p" or "d".
#' @param type type of basis correlation can be legendre or cosine.
#' @param ngrid number of Gauss-Legendre nodes per axis used to evaluate the
#'   double integral. Defaults to 40 for `method = "p"` (accurate to within
#'   about 6e-5 for degree up to `(7, 8)`, validated across several copula
#'   families) and 80 for `method = "d"` (the density integrand converges
#'   more slowly for strongly dependent, low-degrees-of-freedom copulas such
#'   as `tCopula`; push `ngrid` toward 150-200 for those if `"p"`-level
#'   accuracy is needed).
#' @param ... other arguments to function.
#'
#' @return value of polynomial rank correlation.
#' @export
#' @import copula
#'
#' @examples
#' basiscorcopula(copula::claytonCopula(2), 2, 2, copobj = TRUE)
basiscorcopula <- function(cop, j = 1L, k = 1L, copobj, method = "p",
                          type = "legendre", ngrid = if (method == "d") 80L else 40L, ...){
  method <- match.arg(method, c("p", "d"))
  type <- match.arg(type, c("legendre", "cosine"))
  dg <- validate_degrees(j, k)
  if (!is.na(dg$shortcut))
    return(dg$shortcut)
  j <- dg$j
  k <- dg$k
  if (copobj){
    if (!methods::is(cop, "parCopula"))
      stop("Function requires a parametric copula object", call. = FALSE)
    if (methods::is(cop, "tCopula") && cop@parameters[cop@param.names == "df"] %% 1 != 0) # pCopula not implemented for non-integer df
      method <- "d"
  }
  result <- basiscorcopula_gl(cop, j, k, copobj, method, type, ngrid, ...)
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

#' Inner Integrand for Basis Correlation
#'
#' @param u vector argument of function.
#' @param v vector argument of function.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param copobj logical parameter for copula object.
#' @param cop an object of class parCopula or a function.
#' @param method method of calculation which can be "p" or "d".
#' @param type type of basis correlation can be legendre or cosine.
#'
#' @return value of inner integrand
#' @keywords internal
basiscor_inner <- function(u, v, j, k, copobj, cop, method, type, ...){
  if ((copobj) & (method == "p"))
    part1 <- pCopula(cbind(u, v), cop)
  else if ((copobj) & (method == "d"))
    part1 <- dCopula(cbind(u, v), cop)
  else
    part1 <- cop(u, v, ...)
  output <- switch(method,
    d = switch(type,
      legendre = sLegendre(u, j) * sLegendre(v, k) * sqrt(2 * j + 1) * sqrt(2 * k + 1) * part1,
      cosine = 2 * cos(j * pi * u) * cos(k * pi * v) * (-1)^(j + k) * part1,
      stop("Unknown type of basis function", call. = FALSE)
    ),
    p = switch(type,
      legendre = sLegendre(u, j, deriv = TRUE) * sLegendre(v, k, deriv = TRUE) * part1 / (j * (j + 1) * k * (k + 1)),
      cosine = 2 * j * k * (pi^2) * sin(j * pi * u) * sin(k * pi * v) * (-1)^(j + k) * part1,
      stop("Unknown type of basis function", call. = FALSE)
    ),
    stop("Unknown method", call. = FALSE)
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
  type <- match.arg(type, c("legendre", "cosine"))
  method <- match.arg(method, c("T1", "T2", "T3", "T4", "T5", "T6"))
  if (! methods::is(data, "matrix"))
    stop("Must supply data matrix", call. = FALSE)
  if (ncol(data) != 2)
    stop("Data matrix should have 2 columns exactly", call. = FALSE)
  dg <- validate_degrees(j, k)
  if (!is.na(dg$shortcut))
    return(dg$shortcut)
  j <- dg$j
  k <- dg$k
  if ((type != "legendre") & (method == "T6"))
    stop("Discrete polynomials only available for Legendre correlations", call. = FALSE)
  R1 <- rank(data[,1])
  R2 <- rank(data[,2])
  n <- length(R1)
  switch(method,
         T1 = mean(basisfunc(R1/(n+1), j, type) * basisfunc(R2/(n+1), k, type)),
         T2 = mean(basisfunc((R1-0.5)/n, j, type) * basisfunc((R2-0.5)/n, k, type)),
         T3 = stats::cor(basisfunc(R1/(n+1), j, type), basisfunc(R2/(n+1), k, type)),
         T4 = stats::cor(basisfunc((R1-0.5)/n, j, type), basisfunc((R2-0.5)/n, k, type)),
         T5 = n * sum((basisintegral(R1/n, j, type) - basisintegral((R1-1)/n, j, type))*
                      (basisintegral(R2/n, k, type) - basisintegral((R2-1)/n, k, type))),
         T6 = sqrt((2*j+1)*(2*k+1))*exp((lff(n-1,j)+lff(n-1,k)-lff(n+j,j)-lff(n+k,k))/2) *
           mean(discreteLegendre(R1,n,j)*discreteLegendre(R2,n,k)),
         stop("Unknown method", call. = FALSE))
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