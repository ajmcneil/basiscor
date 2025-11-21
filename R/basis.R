# FUNCTIONS FOR GENERAL BASIS FUNCTIONS

#' Orthonormal Basis Function
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

#' Integrated Orthonormal Basis Function
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of integrated basis function.
#' @export
#'
basisintegral <- function(x, degree, type ="legendre"){
  switch(type,
         legendre = ((basisfunc(x, degree + 1, "legendre")/sqrt(2*degree + 3)) -
                       (basisfunc(x, degree - 1, "legendre")/sqrt(2*degree -1)))/(2*sqrt(2*degree + 1)),
         cosine = sin(degree * pi * x) * sqrt(2) * ((-1)^degree)/(degree * pi))
}

#' Derivative of Orthonormal Basis Function
#'
#' @param x vector of values at which basis function should be evaluated.
#' @param degree non-negative integer giving degree of function.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#'
#' @returns vector of values of derivative of basis function.
#' @export
#'
basisderiv <- function(x, degree, type ="legendre"){
  switch(type,
         legendre = sLegendre(x, degree, deriv = TRUE) * sqrt(2 * degree + 1),
         cosine = -degree * pi * sin(degree * pi * x) * sqrt(2) * ((-1)^degree))
}

#' Uniform Distribution Preserving Transformation for Basis Function
#'
#' @param u vector of values between 0 and 1.
#' @param degree degree of polynomial.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#' @param roottol tolerance passed to distLegendre.
#'
#' @return values of udp transformation
#' @export
#'
#' @examples
#' udpfunc(0.3, 3)
udpfunc <- function(u, degree, type ="legendre", roottol = 1e-10){
  switch(type,
         legendre = {
           Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
           distLegendre(stats::predict(Ppoly, u), degree, roottol)$df
           },
         cosine = 1 - acos(cos(degree * pi *u) * (-1)^degree)/pi
         )
}

#' Derivative of Uniform Distribution Preserving Transformation for Basis Function
#'
#' @param u vector of values between 0 and 1.
#' @param degree degree of polynomial.
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#' @param roottol tolerance passed to distLegendre.
#'
#' @return values of derivative of udp transformation
#' @export
#'
#' @examples
#' udpderiv(0.3, 3)
udpderiv <- function(u, degree, type ="legendre", roottol = 1e-10){
  switch(type,
         legendre = {
           Ppoly <- orthopolynom::slegendre.polynomials(degree + 1)[[degree + 1]]
           PpolyD <- stats::deriv(Ppoly, "x")
           distLegendre(stats::predict(Ppoly, u), degree, roottol)$density * stats::predict(PpolyD, u)
         },
         cosine = degree * (-1)^(degree + 1) * sin(degree * pi *u) / abs(sin(degree * pi * u))
  )
}

#' Stochastic Inverse for udp Transformation Associated with Basis Function
#'
#' @param v vector of values between 0 and 1
#' @param degree degree of function
#' @param type character string specifying type of basis function and taking values "legendre" or "cosine".
#' @param roottol tolerance for roots of shifted Legendre polynomial
#' @param endtol tolerance at endpoints for qLegendre
#'
#' @return list consisting of two components - 
#' a list of possible inverse values and a list of probabilities for each value
#' @export
#'
#' @examples
#' udpstochinv(c(0.4,0.5), degree = 4)
udpstochinv <- function(v, degree, type = "legendre", roottol = 1e-10, endtol = 1e-10){
  if (degree == 1)
    return(list(values = as.list(v), probs = as.list(rep(1, length(v)))))
  values <- vector("list", length(v))
  probabilities <- vector("list", length(v))
  if (type == "legendre"){
    y <- qLegendre(v, degree, roottol, endtol)
    legdata <- legendre_initial(degree)
  }
  if (type == "cosine")
    y <- sqrt(2) * cos(pi * (1-v))
  for (i in 1:length(v)){
    rroots <- switch(type,
                     cosine = cosine_roots(y[i], degree),
                     legendre = legendre_roots(legdata$cfs, y[i], roottol))
    values[[i]] <- rroots
    probs <- 1 / abs(udpderiv(rroots, degree, type, roottol))
    probs <- probs/sum(probs) # corrects Legendre case if necessary
    if (is.na(sum(probs))) # in case derivative fails
      probs <- rep(1 / length(probs), length(probs))
    probabilities[[i]] <- probs
  }
  list(values = values, probs = probabilities)
}

#' Find Roots of Cosine Basis Function
#'
#' @param x value to find roots for
#' @param degree degree of function
#'
#' @returns vector of roots
#' @keywords internal
cosine_roots <- function(x, degree = 1){
  theta <- acos(x * (-1)^degree / sqrt(2))
  sign <- ((-1)^(1:degree)) * (-1)
  top <- floor(degree/2)
  add <- c(0, rep(1:top, each = 2))[1:degree]
  (theta*sign + add*2*pi)/(degree*pi)
}

#' Sample From Copula Attaining Basis Correlation Bounds
#'
#' @param j degree of first polynomial
#' @param k degree of second polynomial
#' @param n number of sample points
#' @param type characterstring giving type of basis correlation
#' @param case character variable which should be max or min 
#' @param roottol tolerance for roots of shifted Legendre polynomial
#' @param endtol tolerance at endpoints for qLegendre
#'
#' @return a matrix with nrows and 2 columns containing simulated values from copula
#' @export
#'
#' @examples
#' rstochinvbounds(4, 4)
rstochinvbounds <- function(j, k, n = 1000, type = "legendre", case = "max", roottol = 1e-10, endtol = 1e-10) {
  Ustar <- stats::runif(n)
  Udata <- udpstochinv(Ustar, j, type, roottol, endtol)
  if (case == "min")
    Vdata <- udpstochinv(1-Ustar, k, type, roottol, endtol)
  else if ((case == "max") & (k != j))
    Vdata <- udpstochinv(Ustar, k, type, roottol, endtol)
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

#' Bivariate Stochastic Inversion Using Correlation Basis
#'
#' @param data bivariate matrix of copula data
#' @param j degree of first polynomial
#' @param k degree of second polynomial
#' @param rand randomizer data
#' @param type character string giving type of basis correlation
#' @param roottol tolerance for roots of shifted Legendre polynomial
#' @param endtol tolerance at endpoints for qLegendre
#'
#' @return a matrix with nrows and 2 columns containing simulated values from copula
#' @export
#'
rstochinv <- function(data, j, k, rand = NA, type = "legendre", roottol = 1e-10, endtol = 1e-10) {
  if (!(is.matrix(data)))
    stop("data must be a matrix")
  if (ncol(data) != 2)
    stop("data must be bivariate matrix")
  n <- nrow(data)
  if (is.na(sum(rand)))
    rand <- cbind(stats::runif(n), stats::runif(n))
  if (nrow(rand) != n)
    stop("Randomizer wrong dimension")
  Udata <- udpstochinv(data[,1], j, type, roottol, endtol)
  Vdata <- udpstochinv(data[,2], k, type, roottol, endtol)
  U <- rep(NA, n)
  V <- rep(NA, n)
  for (i in 1:n){
      rroots <- Udata$values[[i]]
      probs <- Udata$probs[[i]]
      nvals <- length(probs)
      cprobs <- c(0, cumsum(probs)[-nvals])
      U[i] <- rroots[max((1:nvals)[rand[i, 1] > cprobs])]
      rroots <- Vdata$values[[i]]
      probs <- Vdata$probs[[i]]
      nvals <- length(probs)
      cprobs <- c(0, cumsum(probs)[-nvals])
      V[i] <- rroots[max((1:nvals)[rand[i, 2] > cprobs])]
  }
  cbind(U, V)
}
