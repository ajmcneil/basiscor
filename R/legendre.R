# FUNCTIONS FOR SHIFTED LEGENDRE POLYNOMIALS

#' Shifted Legendre Polynomial and Derivative
#'
#' @param x vector of values at which polynomial should be evaluated.
#' @param degree non-negative integer giving degree of polynomial.
#' @param deriv logical flag for comuting derivative
#'
#' @return vector of values of polynomial.
#' @export
#'
#' @examples
#' sLegendre(0.5, 2)
sLegendre <- function(x, degree = 1, deriv = FALSE){
  poly <- orthopolynom::legendre.polynomials(degree+1)[[degree+1]]
  if (deriv)
  {
    polyD <- stats::deriv(poly, "x")
    return(2 * stats::predict(polyD, 2*x - 1))
  }
  else
    return(stats::predict(poly, 2*x - 1))
}

#' Collect Data for Distribution of Shifted Legendre Polynomial
#'
#' @param degree of polynomial
#'
#' @returns list with coefficients of polynomial, coefficients of derivative and lower bound of polynomial
#' @keywords internal
legendre_initial <- function(degree){
  Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
  PpolyD <- stats::deriv(Ppoly, "x")
  cfs <- coef(Ppoly)
  cfsD <- coef(PpolyD)
  lbound <- -1
  if (degree %% 2 == 0) {
    tps <- Re(polyroot(cfsD))
    lbound <- min(stats::predict(Ppoly, tps))
  }
  return(list(cfs = cfs, cfsD = cfsD, lbound = lbound))
}

#' Find Real Roots of Shifted Legendre Polynomial
#'
#' @param legcfs coefficients of shifted Legendre polynomial
#' @param x real value in range of polynomial
#' @param roottol tolerance for identification of real roots
#'
#' @returns vector of real roots
#' @keywords internal
legendre_roots <- function(legcfs, x, roottol){
  newcfs <- legcfs
  newcfs[1] <- newcfs[1] - x
  roots <- polyroot(newcfs)
  sort(Re(roots)[abs(Im(roots)) < roottol])
}

#' Distribution of Shifted Legendre Polynomial Transformation of Uniform
#'
#' @param x vector of values 
#' @param degree degree of polynomial
#' @param roottol tolerance for assigning real roots
#'
#' @return list with two vector components containing values of cumulative distribution function and values of density
#' @export
#'
#' @examples
#' distLegendre(0.7, degree = 4)
distLegendre <- function(x, degree, roottol = 1e-10) {
  foutput <- rep(0, length(x))
  Foutput <- foutput
  Foutput[x > 1] <- 1
  legdata <- legendre_initial(degree)
  if (degree == 1){
    foutput[(x >= legdata$lbound) & (x <= 1)] <- 0.5
    Foutput <- pmin(pmax((x+1)/2, 0), 1)
  }
  if (degree == 2){
    foutput[(x >= legdata$lbound) & (x <= 1)] <- 1/((2*x +1)*3)[(x >= legdata$lbound) & (x <= 1)]
    foutput <- sqrt(foutput)
    Foutput <- pmin(pmax(sqrt((2*x + 1)/3), 0) , 1)
  }
  if (degree > 2){
  for (i in 1:length(x)) {
     if (x[i] >= legdata$lbound) {
      rroots <- legendre_roots(legdata$cfs, x[i], roottol)
      nroots <- length(rroots)
      if ((degree %% 2 != 0) & (nroots %% 2 == 0))
        stop("Odd degree and even number of roots")
      if ((degree %% 2 == 0) & (nroots %% 2 != 0))
        stop("Even degree and odd number of roots")
      rrootsD <- 1 / copula::polynEval(legdata$cfsD, rroots)
      if (degree %% 2 != 0){
        rroots <- c(0, rroots)
        rrootsD <- c(0, rrootsD)
      }
      Foutput[i] <- sum(rroots * (-1)^(1:length(rroots)))
      foutput[i] <- sum(rrootsD * (-1)^(1:length(rrootsD)))
    }
  }
  }
  list(df = Foutput, density = foutput)
}

#' Quantile Function of Shifted Legendre Polynomial Transformation of Uniform
#'
#' @param y vector of values 
#' @param degree degree of polynomial
#' @param roottol tolerance for roots of shifted Legendre polynomial
#' @param endtol tolerance at endpoints
#'
#' @return values of quantiles
#' @export
#'
#' @examples
#' qLegendre(0.7, degree = 4)
qLegendre <- function(y, degree, roottol = 1e-10, endtol = 1e-10){
  if (degree == 1)
    return(pmin(pmax((2*y - 1), -1), 1))
  if (degree == 2)
    return(pmin(pmax((3*y^2 - 1)/2, -0.5), 1))
  legdata <- legendre_initial(degree)
  rootfunc <- function(x, y, cfs, degree, roottol){
    rroots <- legendre_roots(cfs, x, roottol)
    nroots <- length(rroots)
    if ((degree %% 2 != 0) & (nroots %% 2 == 0))
      stop("Odd degree and even number of roots")
    if ((degree %% 2 == 0) & (nroots %% 2 != 0))
      stop("Even degree and odd number of roots")
    if (degree %% 2 != 0)
      rroots <- c(0, rroots)
    sum(rroots * (-1)^(1:length(rroots))) - y
  }
  output <- rep(NA,length(y))
  for (i in 1:length(y)){
    if (y[i] <= endtol)
      output[i] <- legdata$lbound
    else if (y[i] >= (1-endtol))
      output[i] <- 1
    else{
      # if (((rootfunc(legdata$lbound, y[i], legdata$cfs, degree, roottol) < 0) & (rootfunc(1, y[i], legdata$cfs, degree, roottol) < 0)) |
      #     ((rootfunc(legdata$lbound, y[i], legdata$cfs, degree, roottol) > 0) & (rootfunc(1, y[i], legdata$cfs, degree, roottol) > 0)))
      #   stop(paste("Root problem at", y[i]))
      output[i] <- stats::uniroot(rootfunc, c(legdata$lbound, 1), y = y[i], cfs = legdata$cfs, 
                                  degree = degree, roottol = roottol)$root
    }
  }
  output
}

#' Calculate Extremal Legendre Correlation
#'
#' @param j degree of first polynomial
#' @param k degree of second polynomial
#' @param case character variable which should be max or min
#'
#' @return value of extremal Legendre correlation
#' @export
#'
#' @examples
#' extremalLegendre(3, 4)
extremalLegendre <- function(j, k, case = "max"){
  integrand <- switch(case,
                      max = function(u, j, k){
                        qLegendre(u, j)*qLegendre(u, k)
                      },
                      min = function(u, j, k){
                        qLegendre(u, j)*qLegendre(1 - u, k)
                      },
                      stop("Unknown case for extremum")
  )
  sqrt(2*j + 1) * sqrt(2*k + 1) * stats::integrate(integrand, 0, 1, j = j, k = k)$value
}