# FUNCTIONS FOR GENERAL BASIS FUNCTIONS
#
# The udp transformation associated with a single degree-j basis function
# (legendre or cosine) -- what used to be basiscor's own udpfunc()/
# udpderiv()/udpstochinv()/cosine_roots()/rstochinv()/rstochinvbounds() -- is
# now udp::udplegendre(j) / udp::udpcosine(j), with udptrans()/udpderiv()/
# udpinverse()/udpsi() as the public interface: e.g. what used to be
# rstochinv(data, j, k) is just
# cbind(udp::udpsi(udp::udplegendre(j), data[, 1]),
#       udp::udpsi(udp::udplegendre(k), data[, 2])).
# basisfunc()/basisderiv()/basisintegral() below evaluate the raw orthonormal
# basis functions themselves, which udp does not expose (it only exposes the
# transformation built on top of them) and which corfunctions.R needs
# directly for the correlation integrands.

#' Shifted Legendre Polynomial and Derivative
#'
#' @param x vector of values at which polynomial should be evaluated.
#' @param degree non-negative integer giving degree of polynomial.
#' @param deriv logical flag for computing derivative
#'
#' @return vector of values of polynomial.
#' @export
#'
#' @examples
#' sLegendre(0.5, 2)
sLegendre <- function(x, degree = 1, deriv = FALSE) {
  degree <- as.integer(degree)
  if (degree < 0) {
    stop("degree must be a non-negative integer")
  }
  t <- 2 * x - 1
  L <- rep(1, length(x))
  Lprev <- rep(0, length(x))
  Ld <- rep(0, length(x))
  Ldprev <- rep(0, length(x))
  if (degree >= 1) {
    for (n in 0:(degree - 1)) {
      Lnext <- ((2 * n + 1) * t * L - n * Lprev) / (n + 1)
      Ldnext <- ((2 * n + 1) * (L + t * Ld) - n * Ldprev) / (n + 1)
      Lprev <- L
      Ldprev <- Ld
      L <- Lnext
      Ld <- Ldnext
    }
  }
  if (deriv) {
    return(2 * Ld)
  }
  L
}

#' Orthonormal Basis Function
#'
#' The degree-`degree` member of an orthonormal basis of `L^2[0, 1]`: the
#' orthonormal shifted Legendre polynomial `P_degree(x) = sqrt(2 degree + 1)
#' L_degree(x)` for `type = "legendre"` (the default), or a cosine basis
#' function for `type = "cosine"`. These are the building blocks
#' \code{\link{basiscor}} and \code{\link{basiscordata}} compute
#' correlations between.
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of basis function.
#' @export
#'
#' @examples
#' basisfunc(seq(0, 1, by = 0.25), 3)
#' basisfunc(seq(0, 1, by = 0.25), 3, type = "cosine")
basisfunc <- function(x, degree, type = "legendre") {
  switch(type,
    legendre = sLegendre(x, degree) * sqrt(2 * degree + 1),
    cosine = cos(degree * pi * x) * sqrt(2) * (-1)^degree
  )
}

#' Integrated Orthonormal Basis Function
#'
#' The indefinite integral (antiderivative vanishing at 0) of
#' \code{\link{basisfunc}}. Used by \code{\link{basiscordata}}'s `"T5"`
#' method, which integrates each basis function exactly over a rank's unit
#' interval rather than evaluating it at a single pseudo-observation.
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of integrated basis function.
#' @export
#'
#' @examples
#' basisintegral(seq(0, 1, by = 0.25), 3)
basisintegral <- function(x, degree, type = "legendre") {
  switch(type,
    legendre = ((basisfunc(x, degree + 1, "legendre") / sqrt(2 * degree + 3)) -
      (basisfunc(x, degree - 1, "legendre") / sqrt(2 * degree - 1))) / (2 * sqrt(2 * degree + 1)),
    cosine = sin(degree * pi * x) * sqrt(2) * ((-1)^degree) / (degree * pi)
  )
}

#' Derivative of Orthonormal Basis Function
#'
#' The derivative of \code{\link{basisfunc}}. Used by
#' \code{\link{basiscorcopula}}'s `method = "p"` integrand.
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of derivative of basis function.
#' @export
#'
#' @examples
#' basisderiv(seq(0, 1, by = 0.25), 3)
basisderiv <- function(x, degree, type = "legendre") {
  switch(type,
    legendre = sLegendre(x, degree, deriv = TRUE) * sqrt(2 * degree + 1),
    cosine = -degree * pi * sin(degree * pi * x) * sqrt(2) * ((-1)^degree)
  )
}
