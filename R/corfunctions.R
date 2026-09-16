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
# integral at a single fixed node count, replacing basiscor_outer()'s nested
# stats::integrate(). Used directly by method = "p" (see basiscorcopula()),
# and as the building block of basiscorcopula_gl_adaptive() for method = "d".
basiscorcopula_gl <- function(cop, j, k, copobj, method, type, ngrid, ...) {
  GL <- gauss_legendre_01(ngrid)
  uv <- expand.grid(u = GL$nodes, v = GL$nodes)
  w <- as.vector(outer(GL$weights, GL$weights))
  vals <- basiscor_inner(uv$u, uv$v, j, k, copobj, cop, method, type, ...)
  sum(w * vals)
}

# method = "p"'s integrand is smooth (built from pCopula and the derivative
# of a fixed-degree polynomial/cosine basis function, never the copula
# density): validated against a 300-node-per-axis reference across several
# copula families at degree up to (7, 8), error <= 6e-5 at a single fixed
# ngrid = 40, comfortably inside stats::integrate()'s own ~1.2e-4 default
# tolerance, at roughly 1.5-2x the speed of the nested adaptive calls it
# replaces.
#
# method = "d" uses the copula density directly, which is far more expensive
# per evaluation for elliptical copulas (normal, t) -- exactly the families
# method = "d" exists for, since their pCopula has no closed form -- so the
# same one-vectorized-batch trick is a much bigger win there when it applies.
# It also fixes a real reliability problem: stats::integrate() on the nested
# density integrand returns NA outright for some entirely ordinary cases
# (e.g. normalCopula(0.3) at j = 6, k = 6), not just extreme ones.
#
# But unlike "p", a single fixed node count is the wrong shape of fix for
# "d": how much resolution the density integrand needs turns out to depend
# on the copula's own dependence structure (how concentrated its density
# is), not predictably on j, k -- a strongly dependent, low-degrees-of-
# freedom t copula is already hard at j = k = 1, while an ordinary one stays
# easy even at high degree. A fixed large ngrid is wasted work on the easy
# majority; a fixed small one is inaccurate on the hard minority.
# basiscorcopula_gl_adaptive() instead starts small (ngrid0, cheap even for
# easy cases -- comparable cost to the old adaptive integrate()) and doubles
# the node count, comparing successive estimates, until they agree to
# reltol or ngridmax is reached, so the cost automatically follows the
# difficulty of the specific (copula, j, k) at hand.
basiscorcopula_gl_adaptive <- function(cop, j, k, copobj, method, type, ngrid0, reltol, ngridmax, ...) {
  n <- ngrid0
  val <- basiscorcopula_gl(cop, j, k, copobj, method, type, n, ...)
  while (n < ngridmax) {
    n2 <- min(2L * n, ngridmax)
    val2 <- basiscorcopula_gl(cop, j, k, copobj, method, type, n2, ...)
    converged <- abs(val2 - val) <= reltol * max(1, abs(val2))
    n <- n2
    val <- val2
    if (converged) {
      break
    }
  }
  val
}

#' Compute Basis Correlation
#'
#' The basis correlation between two degree-`j`, `k` orthonormal basis
#' functions (shifted Legendre polynomials by default, or a cosine basis
#' with `type = "cosine"`) of the two margins of a bivariate distribution.
#' At `j = k = 1` this is exactly Spearman's rho: the orthonormal degree-1
#' shifted Legendre polynomial is \eqn{P_1(u) = \sqrt{3}(2u - 1)}, so
#' \code{basiscor(cop, 1, 1) = 3 E[(2U-1)(2V-1)] = 12 Cov(U, V)}, the usual
#' definition. Higher-degree basis correlations pick up dependence that
#' ordinary correlation and Spearman's rho miss, including non-monotonic
#' dependence.
#'
#' `basiscor()` is a single dispatching entry point: it computes the
#' population value from a parametric copula object or a bivariate copula
#' function (via \code{\link{basiscorcopula}}), or the sample value from a
#' bivariate data matrix (via \code{\link{basiscordata}}).
#'
#' @param object an object of class copula or Copula (from the \pkg{copula}
#'   package), a bicop_dist object (from \pkg{rvinecopulib}), a data matrix,
#'   or a bivariate function describing a copula.
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
#' basiscor(copula::claytonCopula(2), 1, 1) # Spearman's rho
#' data <- copula::rCopula(1000, copula::claytonCopula(2))
#' basiscor(data, 2, 2) # sample analogue, via basiscordata()
#' basiscor(rvinecopulib::bicop_dist("clayton", 0, 2), 2, 2)
basiscor <- function(object, j = 1L, k = 1L, ...){
  if (methods::is(object, "matrix"))
    basiscordata(object, j, k, ...)
  else if ((methods::is(object, "copula")) | (methods::is(object, "Copula")))
    basiscorcopula(object, j, k, copobj = TRUE, ...)
  else if (methods::is(object, "bicop_dist"))
    basiscorcopula(object, j, k, copobj = "bicop", ...)
  else if (methods::is(object, "function"))
    basiscorcopula(object, j, k, copobj = FALSE, ...)
  else
    stop("Unknown object type")
}

#' Compute Basis Correlation for Copula Object or Function
#'
#' The population basis correlation \code{\link{basiscor}} computes from a
#' copula, evaluated by numerical integration of a double integral over the
#' unit square rather than sampling. \code{cop} may be a `parCopula` object
#' from the \pkg{copula} package (\code{copobj = TRUE}), a `bicop_dist`
#' object from the \pkg{rvinecopulib} package (\code{copobj = "bicop"}), or a
#' bivariate function \code{cop(u, v, ...)} giving the copula's own
#' distribution function directly (\code{copobj = FALSE}); the function form
#' is only usable with \code{method = "p"}, since there is no density to
#' fall back on.
#'
#' @param cop an object of class parCopula, an object of class bicop_dist, or a function.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param copobj indicates the kind of copula object `cop` is: `TRUE` for a
#'   `parCopula` object from the \pkg{copula} package, `"bicop"` for a
#'   `bicop_dist` object from \pkg{rvinecopulib}, `FALSE` for a raw bivariate
#'   function.
#' @param method method of calculation: `"p"` integrates against the
#'   copula's distribution function (`copula::pCopula()` or
#'   `rvinecopulib::pbicop()`; the default -- faster, and the only option
#'   when a family has no closed-form density); `"d"` integrates against the
#'   density (`copula::dCopula()` or `rvinecopulib::dbicop()`; needed when a
#'   family has no closed-form distribution function, such as the `copula`
#'   package's `t` copula with non-integer degrees of freedom, which is
#'   detected automatically and switches `method` to `"d"` regardless of
#'   what was requested -- `rvinecopulib`'s `bicop_dist` `t` copula has no
#'   such restriction, so no such switch happens for `copobj = "bicop"`).
#' @param type type of basis correlation can be legendre or cosine.
#' @param ngrid number of Gauss-Legendre nodes per axis used to evaluate the
#'   double integral for `method = "p"`. The default of 40 is accurate to
#'   within about 6e-5 for degree up to `(7, 8)`, validated across several
#'   copula families. Ignored for `method = "d"`; see `ngrid0`.
#' @param ngrid0,reltol,ngridmax for `method = "d"` only: the density
#'   integrand can need anywhere from very few to very many nodes depending
#'   on the copula's own dependence structure, not predictably on `j`, `k`,
#'   so the node count is chosen adaptively -- starting at `ngrid0` nodes per
#'   axis, doubling and comparing successive estimates until they agree to
#'   within `reltol` (relative to the finer one) or `ngridmax` is reached.
#'   Defaults `ngrid0 = 20`, `reltol = 1e-4`, `ngridmax = 320`.
#' @param ... other arguments to function.
#'
#' @return value of polynomial rank correlation.
#' @export
#' @import copula
#' @import rvinecopulib
#'
#' @examples
#' basiscorcopula(copula::claytonCopula(2), 2, 2, copobj = TRUE)
#' basiscorcopula(rvinecopulib::bicop_dist("clayton", 0, 2), 2, 2, copobj = "bicop")
basiscorcopula <- function(cop, j = 1L, k = 1L, copobj, method = "p",
                          type = "legendre", ngrid = 40L,
                          ngrid0 = 20L, reltol = 1e-4, ngridmax = 320L, ...){
  method <- match.arg(method, c("p", "d"))
  type <- match.arg(type, c("legendre", "cosine"))
  dg <- validate_degrees(j, k)
  if (!is.na(dg$shortcut))
    return(dg$shortcut)
  j <- dg$j
  k <- dg$k
  if (isTRUE(copobj)){
    if (!methods::is(cop, "parCopula"))
      stop("Function requires a parametric copula object", call. = FALSE)
    if (methods::is(cop, "tCopula") && cop@parameters[cop@param.names == "df"] %% 1 != 0) # pCopula not implemented for non-integer df
      method <- "d"
  } else if (identical(copobj, "bicop")){
    if (!methods::is(cop, "bicop_dist"))
      stop("Function requires a bicop_dist object from rvinecopulib", call. = FALSE)
  }
  result <- if (method == "p") {
    basiscorcopula_gl(cop, j, k, copobj, method, type, ngrid, ...)
  } else {
    basiscorcopula_gl_adaptive(cop, j, k, copobj, method, type, ngrid0, reltol, ngridmax, ...)
  }
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
#' @param copobj indicates the kind of copula object `cop` is: `TRUE`, `"bicop"` or `FALSE` (see \code{\link{basiscorcopula}}).
#' @param cop an object of class parCopula, an object of class bicop_dist, or a function.
#' @param method method of calculation which can be "p" or "d".
#' @param type type of basis correlation can be legendre or cosine.
#'
#' @return value of inner integrand
#' @keywords internal
basiscor_inner <- function(u, v, j, k, copobj, cop, method, type, ...){
  if (isTRUE(copobj) && method == "p")
    part1 <- pCopula(cbind(u, v), cop)
  else if (isTRUE(copobj) && method == "d")
    part1 <- dCopula(cbind(u, v), cop)
  else if (identical(copobj, "bicop") && method == "p")
    part1 <- pbicop(cbind(u, v), cop)
  else if (identical(copobj, "bicop") && method == "d")
    part1 <- dbicop(cbind(u, v), cop)
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
#' Calls \code{\link{basiscor}} for every pair of degrees `1:maxorder`,
#' returning the results as a matrix -- entry `[j, k]` is `basiscor(object,
#' j, k, ...)`. This is what \code{\link{basisexpand}} uses to find its
#' basis expansion.
#'
#' @param object an object of class copula or Copula, or a data matrix, or a bivariate function describing a copula.
#' @param maxorder maximum order of the polynomials.
#' @param symmetric logical variable stating whether the matrix is known a
#'   priori to be symmetric (true, for instance, for an exchangeable copula
#'   with `j`, `k` given by the same margin ordering both times), halving
#'   the number of calls to \code{\link{basiscor}}.
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
#' The sample analogue of \code{\link{basiscor}}: a generalization of
#' Spearman's rho computed from the ranks of a bivariate sample. At `j = k =
#' 1` with the default `method = "T3"`, it is exactly the ordinary sample
#' Spearman correlation.
#'
#' @param data a matrix of data wit two columns.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#' @param method method of calculation, turning the ranks `R1`, `R2` of the
#'   two columns (sample size `n`) into a single value:
#'   \describe{
#'     \item{`"T1"`}{`mean(basisfunc(R1/(n+1), j) * basisfunc(R2/(n+1), k))`,
#'       a covariance-type statistic on the midpoint pseudo-observations
#'       `R/(n+1)`.}
#'     \item{`"T2"`}{as `"T1"`, but on pseudo-observations `(R-0.5)/n`.}
#'     \item{`"T3"`}{the default; `stats::cor()` (rather than the raw mean
#'       product) of the same basis-function values as `"T1"`.}
#'     \item{`"T4"`}{`stats::cor()` of the same values as `"T2"`.}
#'     \item{`"T5"`}{exact integral of the basis functions over each rank's
#'       unit interval, rather than evaluating at a single pseudo-observation.}
#'     \item{`"T6"`}{an exact discrete-Legendre-polynomial estimator;
#'       `type = "legendre"` only.}
#'   }
#'
#' @return sample polynomial rank correlation value.
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::claytonCopula(2))
#' basiscordata(data, 1, 1) # equals cor(data[, 1], data[, 2], method = "spearman")
#' basiscordata(data, 2, 3)
#'
basiscordata <- function(data, j, k, type = "legendre", method = "T3"){
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
#' @param n size of sample.
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