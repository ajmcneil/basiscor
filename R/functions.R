#' Shifted Legendre Polynomial
#'
#' @param x vector of values at which polynomial should be evaluated.
#' @param degree non-negative integer giving degree of polynomial.
#'
#' @return vector of values of polynomial.
#' @export
#'
#' @examples
#' sLegendre(0.5, 2)
sLegendre <- function(x, degree = 1){
  poly <- orthopolynom::legendre.polynomials(degree+1)[[degree+1]]
  stats::predict(poly, 2*x - 1)
}

#' Orthonormal basis functions
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of basis function.
#' @export
#'
basisfunc <- function(x, degree, type ="legendre"){
  switch(type,
         legendre = sLegendre(x, degree) * sqrt(2 * degree + 1),
         cosine = cos(degree * pi * x) * sqrt(2) * (-1)^degree)
}

#' Derivative of Shifted Legendre Polynomial
#' 
#' The function with non-negative degree n is the derivative of the shifted Legendre of degree n + 1.
#'
#' @param x vector of values at which polynomial should be evaluated.
#' @param degree non-negative integer giving degree of polynomial.
#'
#' @return vector of values of polynomial.
#' @export
#'
#' @examples
#' sLegendreD(0.5, 2)
sLegendreD <- function(x, degree = 1){
  poly <- orthopolynom::legendre.polynomials(degree+1)[[degree+1]]
  polyD <- stats::deriv(poly, "x")
  2 * stats::predict(polyD, 2*x - 1)
}

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
  if (type == "legendre")
    output <- switch(method,
                     "p" = sqrt(2 * j + 1) * sqrt(2 * k + 1) * (result*j*(j+1)*k*(k+1) - 1),
                     "d" = sqrt(2 * j + 1) * sqrt(2 * k + 1) * result)
  if (type == "cosine")
    output <- switch(method,
                     "p" = ((-1)^(j+k))* 2 * result * j * k * pi^2 - 2,
                     "d" =  ((-1)^(j+k))* 2 * result)
  output
}

#' Outer integrand for basis correlation
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

#' Inner integrand for basis correlation
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
  if (type == "cosine")
    output <- switch(method,
         "p" = sin(j*pi*u) * sin(k*pi*v) * part1,
         "d" = cos(j*pi*u) * cos(k*pi*v) * part1)
  if (type == "legendre")
    output <- switch(method,
  "p" = sLegendreD(u,j) * sLegendreD(v, k) * part1 / (j*(j+1)*k*(k+1)),
  "d" = sLegendre(u,j) * sLegendre(v, k) * part1)
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
#' basiscorM(copula::claytonCopula(2))
basiscorM <- function(object, maxorder = 4, symmetric = FALSE, ...){
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



#' Distribution function of Legendre polynomial transformation of uniform
#'
#' @param x vector of values 
#' @param degree degree of polynomial
#' @param tol tolerance for assigning real roots
#'
#' @return values of cumulative distribution function
#' @export
#'
#' @examples
#' pLegendre(0.7, degree = 4)
pLegendre <- function(x, degree, tol = 1e-10) {
  if (degree == 1)
    return(pmin(pmax((x+1)/2, 0), 1))
  if (degree == 2)
    return(pmin(pmax(sqrt((2*x + 1)/3), 0) , 1))
  Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
  PpolyD <- stats::deriv(Ppoly, "x")
  cfs <- coef(Ppoly)
  cfsD <- coef(PpolyD)
  degreeodd <-  (degree %% 2 != 0)
  lbound <- -1
  if (!degreeodd) {
    tps <- Re(polyroot(cfsD))
    lbound <- min(stats::predict(Ppoly, tps))
  }
  Foutput <- rep(0, length(x))
  for (i in 1:length(x)) {
    if (x[i] > 1) 
      Foutput[i] <- 1
    else if (x[i] >= lbound) {
      newcfs <- cfs
      newcfs[1] <- newcfs[1] - x[i]
      roots <- polyroot(newcfs)
      rroots <- sort(Re(roots)[abs(Im(roots)) < tol])
      nroots <- length(rroots)
      if (degreeodd & (nroots %% 2 == 0))
        stop("Odd degree and even number of roots")
      if ((!degreeodd) & (nroots %% 2 != 0))
        stop("Even degree and odd number of roots")
      if (degreeodd) 
        rroots <- c(0, rroots)
      Foutput[i] <- sum(rroots * (-1)^(1:length(rroots)))
    }
  }
  Foutput
}

#' Density of Legendre polynomial transformation of uniform
#'
#' @param x vector of values 
#' @param degree degree of polynomial
#' @param tol tolerance for assigning real roots
#'
#' @return values of density
#' @export
#'
#' @examples
#' dLegendre(0.7, degree = 4)
dLegendre <- function(x, degree, tol = 1e-10) {
  Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
  PpolyD <- stats::deriv(Ppoly, "x")
  cfs <- coef(Ppoly)
  cfsD <- coef(PpolyD)
  degreeodd <-  (degree %% 2 != 0)
  lbound <- -1
  if (!degreeodd) {
    tps <- Re(polyroot(cfsD))
    lbound <- min(stats::predict(Ppoly, tps))
  }
  foutput <- rep(0, length(x))
  if (degree == 1)
    foutput[(x >= lbound) & (x <= 1)] <- 0.5
  if (degree == 2){
    foutput[(x >= lbound) & (x <= 1)] <- 1/((2*x +1)*3)[(x >= lbound) & (x <= 1)]
    foutput <- sqrt(foutput)
  }
  if (degree > 2){
    for (i in 1:length(x)) {
      if ((x[i] >= lbound) & (x[i] <= 1)) {
        newcfs <- cfs
        newcfs[1] <- newcfs[1] - x[i]
        roots <- polyroot(newcfs)
        rroots <- sort(Re(roots)[abs(Im(roots)) < tol])
        nroots <- length(rroots)
        if (degreeodd & (nroots %% 2 == 0))
          stop("Odd degree and even number of roots")
        if ((!degreeodd) & (nroots %% 2 != 0))
          stop("Even degree and odd number of roots")
        rrootsD <- 1 / stats::predict(PpolyD, rroots)
        if (degreeodd) 
          rrootsD <- c(0, rrootsD)
        foutput[i] <- sum(rrootsD * (-1)^(1:length(rrootsD)))
      }
    }
  }
  foutput
}

#' Quantile function of Legendre polynomial transformation of uniform
#'
#' @param y vector of values 
#' @param degree degree of polynomial
#' @param tol tolerance at endpoints
#'
#' @return values of quantiles
#' @export
#'
#' @examples
#' qLegendre(0.7, degree = 4)
qLegendre <- function(y, degree, tol = 1e-10){
  if (degree == 1)
    return(pmin(pmax((2*y - 1), -1), 1))
  if (degree == 2)
    return(pmin(pmax((3*y^2 - 1)/2, -0.5), 1))
  Ppoly <- orthopolynom::slegendre.polynomials(degree+1)[[degree+1]]
  PpolyD <- stats::deriv(Ppoly, "x")
  lbound <- -1
  if (degree %% 2 == 0){
    tps <- Re(polyroot(coef(PpolyD)))
    lbound <- min(stats::predict(Ppoly, tps))
  }
  rootfunc <- function(x, y, degree){
    pLegendre(x, degree) - y
  }
  output <- rep(NA,length(y))
  for (i in 1:length(y)){
    if (y[i] <= tol)
      output[i] <- lbound
    else if (y[i] >= (1-tol))
      output[i] <- 1
    else{
      if (((rootfunc(lbound, y[i], degree) < 0) & (rootfunc(1, y[i], degree) < 0)) |
          ((rootfunc(lbound, y[i], degree) > 0) & (rootfunc(1, y[i], degree) > 0)))
        stop(paste("Root problem at", y[i]))
      output[i] <- stats::uniroot(rootfunc, c(lbound, 1), y = y[i],degree = degree)$root
    }
  }
  output
}

#' Uniform distribution preserving transformation for Legendre polynomial
#'
#' @param u vector of values between 0 and 1
#' @param degree degree of polynomial
#' @param tol tolerance passed to pLegendre
#'
#' @return values of udp transformation
#' @export
#'
#' @examples
#' udpLegendre(0.3, 3)
udpLegendre <- function(u, degree, tol = 1e-10) {
  Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
  x <- stats::predict(Ppoly, u)
  pLegendre(x, degree, tol = 1e-10)
}

#' Density of udp transformation for Legendre polynomial
#'
#' @param u vector of values between 0 and 1
#' @param degree degree of polynomial
#' @param tol tolerance passed to dLegendre
#'
#' @return values of density of udp transformation
#' @export
#'
#' @examples
#' udpLegendreD(0.3, 3)
udpLegendreD <- function(u, degree, tol = 1e-10) {
  Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
  PpolyD <- stats::deriv(Ppoly, "x")
  x <- stats::predict(Ppoly, u)
  dLegendre(x, degree, tol = 1e-10) * stats::predict(PpolyD, u)
}

#' Data for stochastic inverse of udp Legendre transformation
#'
#' @param v vector of values between 0 and 1
#' @param degree degree of polynomial
#' @param tol tolerance passed to dLegendre and used in root finding
#' @param data logical value specifying whether data for stochastic inverse should be returned
#' @param Z the randomizer uniform variable which should be same length as v
#'
#' @return either the inverse values or a list consisting of two components - 
#' a list of possible inverse values and a list of probabilities for each value
#' @export
#'
#' @examples
#' udpLegendre(c(0.4,0.5), degree = 4)
udpLegendreI <- function(v, degree, tol = 1e-10, data = FALSE, 
                         Z = stats::runif(length(v))) {
if (degree == 1)
  return(list(values = as.list(v), probs = as.list(rep(1, length(v)))))
Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
cfs <- coef(Ppoly)
lbound <- -1
if (degree %% 2 == 0){
  tps <- Re(polyroot(coef(stats::deriv(Ppoly, "x"))))
  lbound <- min(stats::predict(Ppoly, tps))
}
rootfunc <- function(x, v, degree) {
  pLegendre(x, degree) - v
}
if (data){
  values <- vector("list", length(v))
  probabilities <- vector("list", length(v))
}
else
  output <- rep(NA, length(v))
for (i in 1:length(v)) {
  if (v[i] <= tol)
    y <- lbound
  else if (v[i] >= (1-tol))
    y <- 1
  else
    y <- stats::uniroot(rootfunc, c(lbound, 1), v = v[i], degree = degree)$root
  newcfs <- cfs
  newcfs[1] <- newcfs[1] - y
  roots <- polyroot(newcfs)
  rroots <- Re(roots)[abs(Im(roots)) < tol]
  probs <- 1 / abs(udpLegendreD(rroots, degree, tol))
  probs <- probs / sum(probs)
  if (is.na(sum(probs)))
    probs <- rep(1 / length(probs), length(probs))
  if (data){
    values[[i]] <- rroots
    probabilities[[i]] <- probs
  }
  else
  {
    nvals <- length(probs)
    cprobs <- c(0,cumsum(probs)[-nvals])
    output[i] <- rroots[max((1:nvals)[Z[i] > cprobs])]
  }
}
if (data)
  return(list(values = values, probs = probabilities))
return(output)
}

#' Sample from copula attaining extremal Legendre correlation
#'
#' @param j degree of first polynomial
#' @param k degree of second polynomial
#' @param n number of sample points
#' @param case character variable which should be max or min 
#'
#' @return a matrix with nrows and 2 columns containing simulated values from copula
#' @export
#'
#' @examples
#' rextremalLegendre(4, 4)
rextremalLegendre <- function(j, k, n = 1000, case = "max") {
  Ustar <- stats::runif(n)
  Udata <- udpLegendreI(Ustar, j, data = TRUE)
  if (case == "min")
    Vdata <- udpLegendreI(1-Ustar, k, data = TRUE)
  else if ((case == "max") & (k != j))
    Vdata <- udpLegendreI(Ustar, k, data = TRUE)
  else 
    Vdata <- Udata
  U <- rep(NA, n)
  V <- rep(NA, n)
  for (i in 1:n){
    choiceU <- 1
    choiceV <- 1
    # multinomial sampling is independent but need not be
    if (j > 1)
      choiceU <- which(as.vector(stats::rmultinom(1, 1, Udata$probs[[i]])) == 1)
    U[i] <- Udata$values[[i]][choiceU]
    if (k >1)
      choiceV <- which(as.vector(stats::rmultinom(1, 1, Vdata$probs[[i]])) == 1)
    V[i] <- Vdata$values[[i]][choiceV]
  }
  cbind(U, V)
}

#' Calculate extremal Legendre correlation
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
         T5 = n * sum((basisfuncI(R1/n, j, type) - basisfuncI((R1-1)/n, j, type))*
                      (basisfuncI(R2/n, k, type) - basisfuncI((R2-1)/n, k, type))),
         T6 = sqrt((2*j+1)*(2*k+1))*exp((lff(n-1,j)+lff(n-1,k)-lff(n+j,j)-lff(n+k,k))/2) *
           mean(discreteLegendre(R1,n,j)*discreteLegendre(R2,n,k)))
}


#' Compute Sample Polynomial Rank Correlation
#'
#' @param data a matrix of data wit two columns.
#' @param j non-negative integer giving order of first polynomial.
#' @param k non-negative integer giving order of second polynomial.
#' @param type method of calculation
#'
#' @return sample polynomial rank correlation value.
#' @export
#'
polycordata <- function(data, j, k, type = "T1"){
  if (! methods::is(data, "matrix"))
    stop("Must supply data matrix")
  if (ncol(data) != 2)
    stop("Data matrix should have 2 columns exactly")
  R1 <- rank(data[,1])
  R2 <- rank(data[,2])
  n <- length(R1)
  switch(type,
         T1 = sqrt((2*j+1)*(2*k+1))*mean(sLegendre(R1/(n+1),j)*sLegendre(R2/(n+1),k)),
         T2 = sqrt((2*j+1)*(2*k+1))*mean(sLegendre((R1-0.5)/n,j)*sLegendre((R2-0.5)/n,k)),
         T3 = stats::cor(sLegendre(R1/(n+1),j),sLegendre(R2/(n+1),k)),
         T4 = stats::cor(sLegendre((R1-0.5)/n,j),sLegendre((R2-0.5)/n,k)),
         T5 = n*sqrt((2*j+1)*(2*k+1))*sum((sLegendreI(R1/n,j)-sLegendreI((R1-1)/n,j))*
                                            (sLegendreI(R2/n,k)-sLegendreI((R2-1)/n,k))),
         T6 = sqrt((2*j+1)*(2*k+1))*exp((lff(n-1,j)+lff(n-1,k)-lff(n+j,j)-lff(n+k,k))/2) *
           mean(discreteLegendre(R1,n,j)*discreteLegendre(R2,n,k)))
}

#' Integral of Shifted Legendre Polynomial
#' 
#'
#' @param x vector of values at which polynomial should be evaluated.
#' @param degree non-negative integer giving degree of polynomial.
#'
#' @return vector of values of integrated polynomial.
#' @export
#'
#' @examples
#' sLegendreI(0.5, 2)
sLegendreI <- function(x, degree){
  (sLegendre(x,degree+1) - sLegendre(x,degree-1))/(2*(2*degree+1))
}

#' Integrated orthonormal basis functions
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of basis function.
#' @export
#'
basisfuncI <- function(x, degree, type ="legendre"){
  switch(type,
         legendre = ((basisfunc(x, degree + 1, "legendre")/sqrt(2*degree + 3)) -
           (basisfunc(x, degree - 1, "legendre")/sqrt(2*degree -1)))/(2*sqrt(2*degree + 1)),
         cosine = sin(degree * pi * x) * sqrt(2) * ((-1)^degree)/(degree * pi))
}

#' Log of falling factorial function
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