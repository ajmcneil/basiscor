#' Find Basis Expansion Maximizing Generalized Spearman Correlation
#'
#' @param data matrix of bivariate data
#' @param maxorder order of square matrix of basis correlations to be used in expansion
#' @param type character variable specifying type of basis correlation taking values "legendre" or "cosine"
#'
#' @returns a bex object comprising 4 components: type of basis correlation, maximum value of 
#' generalized Spearman correlation, weights for first function, and weights for second function.
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::tCopula(0.5, df = 2))
#' basisexpand(data, maxorder = 8, type = "legendre")
basisexpand <- function(data, maxorder = 8, type= "legendre"){
  basiscormat <- basiscormatrix(data, maxorder = maxorder, type = type)
  svdres <- svd(basiscormat, 1, 1)
  maxcor <- svdres$d[1]
  alphag <- svdres$u
  alphah <- svdres$v
  pcsm <- pcsm_analysis(alphag, type, plotfunc = FALSE)
  if (pcsm$Rpartition[pcsm$M] > pcsm$Rpartition[pcsm$M + 1]){# test for decreasing in final interval of partition
    alphag <- -alphag
    alphah <- -alphah
  }
  results <- list(type=type,
                  maxcor = maxcor, 
                  alphag = alphag,
                  alphah = alphah)
  class(results) <- c("bex", class(results))
  return(results)
}

#' Evaluate Basis Expansion Function
#'
#' @param u vector of values at which to compute expansion
#' @param expansion a basis expansion given by basisexpand
#' @param margin number (1 or 2) specifying the first of second basis expansion function
#'
#' @returns vector of values of the function
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::tCopula(0, df = 1))
#' bes <- basisexpand(data, maxorder = 8, type = "legendre")
#' bexfunc(c( 0.2,0.4,0.6), bes, 1)
bexfunc <- function(u, expansion, margin){
  alpha <- switch(margin,
                  expansion$alphag,
                  expansion$alphah,
                  stop("No such margin"))
  bexfunc_uni(u, alpha, expansion$type)
}

#' Evaluate Distribution Function of Basis Expansion of Uniform
#'
#' @param x vector of values at which to compute distribution function
#' @param expansion a basis expansion given by basisexpand
#' @param margin number (1 or 2) specifying the first or second basis expansion function
#' @param n_grid number of grid points for finding turning points and roots of basis expansion function
#' @param tol tolerance for finding roots
#' @param h tolernace for finding numerical derivatives
#'
#' @returns vector of values of the function
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::tCopula(0, df = 1))
#' bes <- basisexpand(data, maxorder = 8, type = "legendre")
#' pbex(c(-0.5, 4), bes, 1)
pbex <- function(x, expansion, margin, n_grid = 500, tol = 1e-8, h = 1e-6){
  alpha <- switch(margin,
                  expansion$alphag,
                  expansion$alphah,
                  stop("No such margin"))
  type <- expansion$type
  pcsm <- pcsm_analysis(alpha, type, n_grid = n_grid, plotfunc = FALSE)
  M <- pcsm$M
  output <- rep(NA, length(x))
  for (i in 1:length(x)){
    rootanalysis <- root_analysis(x[i], alpha, type,
                                  n_grid = n_grid, tol = tol, h = h)
    output[i] <- sum(rootanalysis$roots * rootanalysis$signs)
    if ((pcsm$Rpartition[M] > x[i]) & (x[i] > pcsm$Rpartition[M+1])) #D_x
      output[i] <- output[i] + 1
    else if (x[i] >= max(pcsm$Rpartition[M], pcsm$Rpartition[M+1])) #O_x
      output[i] <- output[i] + 1
  }
  output
}

#' Evaluate Udp Function for Basis Expansion
#'
#' @param u vector of values at which to compute udp function
#' @param expansion a basis expansion given by basisexpand
#' @param margin number (1 or 2) specifying the first or second basis expansion function
#' @param n_grid number of grid points for finding turning points and roots of basis expansion function
#' @param tol tolerance for finding roots
#' @param h tolerance for finding numerical derivatives
#'
#' @returns vector of values of the function
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::tCopula(0, df = 1))
#' bes <- basisexpand(data, maxorder = 8, type = "legendre")
#' udpbex(c(0.2, 0.4, 0.6), bes, 1)
udpbex <- function(u, expansion, margin, n_grid = 500, tol = 1e-8, h = 1e-6){
  x <- bexfunc(u, expansion, margin)
  pbex(x, expansion, margin, n_grid, tol, h)
}

#' Transform Data with Udp Function for Basis Expansion
#'
#' @param Udata matrix of data with values in (0,1)
#' @param expansion a basis expansion given by basisexpand
#' @param boundaryadjust logical variable specifying whether to keep transformed data strictly
#' in (0,1)
#' @param btol user-specified adjustment for boundaryadjust
#' @param n_grid number of grid points for finding turning points and roots of basis expansion function
#' @param tol tolerance for finding roots
#' @param h tolerance for finding numerical derivatives
#'
#' @returns vector of values of the function
#' @export
udpbexdata <- function(Udata, expansion, boundaryadjust = TRUE, 
                       btol = NA, n_grid = 500, tol = 1e-8, h = 1e-6){
  Vdata1 <- udpbex(Udata[,1], expansion, 1, n_grid, tol, h)
  Vdata2 <- udpbex(Udata[,2], expansion, 2, n_grid, tol, h)
  if (boundaryadjust){
    if (is.na(btol))
      btol <- 1/(2*length(Vdata1))
    Vdata1[Vdata1 == 0] = btol
    Vdata1[Vdata1 == 1] = 1 - btol
    Vdata2[Vdata2 == 0] = btol
    Vdata2[Vdata2 == 1] = 1 - btol
  }
  cbind(Vdata1, Vdata2)
}

#' Plot Basis Expansion
#'
#' @method plot bex
#' @param x a basis expansion object created by basisexpand
#' @param n_grid number of grid points for plotting function
#' @param udp logical value specifying whether to plot udp functions rather than actual basis expansion
#' @param ... other parameters
#' functions
#' @export
plot.bex <- function(x, n_grid = 200, udp = FALSE, ...){
  uvals <- seq(from = 0, to = 1, length = n_grid)
  if (udp){
    gvals <- udpbex(uvals, x, 1)
    hvals <- udpbex(uvals, x, 2)
    yl <- "T(u)"
  }
  else{
    gvals <- bexfunc(uvals, x, 1)
    hvals <- bexfunc(uvals, x, 2)
    yl <- "g(u)"
  }
  plot(uvals, gvals, type = "l", ylab = yl, xlab = "u", ylim = range(gvals,hvals))
  graphics::lines(uvals, hvals, col = "red")
}

#' Evaluate Single Basis Expansion Function
#'
#' @param u vector of values at which to compute expansion
#' @param alpha vector of weights
#' @param type character string specifying type of basis functions, "legendre" or "cosine"
#' @param k additional argument for root finding
#'
#' @keywords internal 
bexfunc_uni <- function(u, alpha, type, k = 0){
  output <- rep(0, length(u))
  for (i in 1:length(alpha))
    output <- output + alpha[i] * basisfunc(u, degree = i, type = type)
  output - k
}

#' Analysis of Piecewise Strictly Monotonic Function
#'
#' @param alpha vector of weights
#' @param type character string specifying type of basis functions, "legendre" or "cosine"
#' @param n_grid number of grid points for finding turning points
#' @param plotfunc logical variable specifying whether a plot should be made
#' 
#' keywords internal
pcsm_analysis <- function(alpha, type, n_grid = 500, plotfunc = FALSE){
  x <- seq(0, 1, length.out = n_grid)
  y <- bexfunc_uni(x, alpha, type)
  xn <- x[-1]
  yn <- sign(diff(y))
  xnn <- xn[-1]
  nsignchanges <- sum(diff(sign(yn)) != 0)
  if (nsignchanges > 0)
    tps <- xnn[which(diff(sign(yn)) != 0)]
  else
    tps <- 0
  if (plotfunc){
    plot(x, y, type="l", xlab="u", ylab="g(u)")
    graphics::abline(v=tps[tps>0], lty= 3, col = "red")
  }
  Apartition <- c(0, tps, 1)
  Rpartition <- bexfunc_uni(Apartition, alpha, type)
  list(M = length(Apartition) - 1, Apartition = Apartition, Rpartition = Rpartition)
}

#' Analysis of Roots of Piecewise Strictly Monotonic Function
#'
#' @param k real scalar value to find roots for
#' @param alpha vector of weights
#' @param type character string specifying type of basis functions, "legendre" or "cosine"
#' @param n_grid number of grid points for finding turning points
#' @param tol tolerance for finding roots
#' @param h tolerance for finding numerical derivatives
#' 
#' keywords internal
root_analysis <- function(k, alpha, type, n_grid = 500, tol = 1e-8, h = 1e-6) {
  # Step 1: Sample the function
  x <- seq(0, 1, length.out = n_grid)
  y <- bexfunc_uni(x, alpha, type, k)
  # Step 2: Detect sign changes
  sign_changes <- which(diff(sign(y)) != 0)
  roots <- numeric(0)
  # Step 3: Locate each root accurately with uniroot()
  for (i in sign_changes) {
    lower <- x[i]
    upper <- x[i + 1]
    if (is.finite(y[i]) && is.finite(y[i + 1]) && y[i] * y[i + 1] <= 0) {
      r <- tryCatch(stats::uniroot(bexfunc_uni, 
                            c(lower, upper), 
                            tol = tol,
                            k = k,
                            alpha = alpha,
                            type = type)$root,
                    error = function(e) NA)
      if (!is.na(r)) roots <- c(roots, r)
    }
  }
  # Step 4: Clean up duplicates (tangential zeros, etc.)
  roots <- sort(unique(round(roots, digits = ceiling(-log10(tol)))))
  roots <- roots[(roots > 0) & (roots < 1)] # must be in (0,1)
  # Step 5: Estimate derivative near each root to determine monotonicity
  signs <- logical(length(roots))
  for (i in seq_along(roots)) {
    xi <- roots[i]
    # central difference derivative
    dfdx <- (bexfunc_uni(xi + h, alpha, type) - 
               bexfunc_uni(xi - h, alpha, type)) / (2 * h)
    signs[i] <- 2*as.numeric(dfdx > 0) - 1
  }
  # Step 6: Return structured output
  list(
    n_roots = length(roots),
    roots = roots,
    signs = signs
  )
}

