# BASIS EXPANSIONS MAXIMIZING GENERALIZED SPEARMAN CORRELATION
#
# The heavy lifting for a linear combination g(u) = sum_i alpha_i *
# basisfunc(u, i, type) -- finding its turning points, the distribution F of
# g(U), the exact pre-images and selection probabilities of T(u) = F(g(u)),
# and drawing the resulting graph -- used to be basiscor's own
# pcsm_analysis()/root_analysis()/bexfunc_uni() (a grid-search-then-uniroot()
# approach). All of that is now udp::udplegendrebex(alpha) /
# udp::udpcosinebex(alpha): an exact, panel-spline construction of the same
# transformation, with udptrans()/udpinverse()/udpsi()/pcoincide()/plot() as
# the public interface. What remains genuinely basiscor's own is choosing
# alpha (an SVD of the basis-correlation matrix) and a sign convention for
# it; the bex class below just packages the two weight vectors together with
# the two udp transformation objects built from them.
#
# basiscor's earlier standalone pbex() -- the distribution function F applied
# to an arbitrary raw value, independent of any u -- has no equivalent public
# entry point in udp (only T = F o g is exposed) and is dropped: any x of
# interest is g(u) for some u in [0, 1], so udp::udptrans(x@Tg, u) already
# covers it.

#' Class of Basis Expansion Objects
#'
#' A `bex` object packages the pair of basis-function expansions `g(u) =
#' sum_i alphag[i] * basisfunc(u, i, type)` and `h(u) = sum_i alphah[i] *
#' basisfunc(u, i, type)` that jointly maximize the generalized Spearman
#' correlation `maxcor` between the two margins of a bivariate sample (see
#' [basisexpand()]), together with the corresponding \pkg{udp}
#' uniform-distribution-preserving transformations `Tg`, `Th` -- `Tg(u) =
#' F(g(u))`, `Th(u) = F(h(u))` -- built from them.
#'
#' @slot type character; `"legendre"` or `"cosine"`.
#' @slot maxcor numeric; the maximal generalized Spearman correlation attained.
#' @slot alphag,alphah numeric; the weights of the two expansions.
#' @slot Tg,Th objects of class \code{udp} (see `udp::udplegendrebex()` /
#'   `udp::udpcosinebex()`) built from `alphag`, `alphah`.
#'
#' @seealso [basisexpand()] to construct one.
#' @importClassesFrom udp udp
#' @export
setClass("bex",
  slots = list(
    type = "character", maxcor = "numeric",
    alphag = "numeric", alphah = "numeric",
    Tg = "udp", Th = "udp"
  )
)

#' @describeIn bex-class Show method for bex objects.
#' @param object an object of class \linkS4class{bex}.
#' @export
setMethod("show", "bex", function(object) {
  cat("An object of class \"bex\"\n")
  cat("type: ", object@type, "  degree: ", length(object@alphag), "\n", sep = "")
  cat("maxcor: ", format(object@maxcor, digits = 4), "\n\n", sep = "")
  cat("alphag:\n")
  print(object@alphag)
  cat("alphah:\n")
  print(object@alphah)
})

#' Find Basis Expansion Maximizing Generalized Spearman Correlation
#'
#' @param data matrix of bivariate data
#' @param maxorder order of square matrix of basis correlations to be used in expansion
#' @param type character variable specifying type of basis correlation taking values "legendre" or "cosine"
#'
#' @returns an object of class \linkS4class{bex}.
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::tCopula(0.5, df = 2))
#' basisexpand(data, maxorder = 8, type = "legendre")
basisexpand <- function(data, maxorder = 8, type = "legendre") {
  basiscormat <- basiscormatrix(data, maxorder = maxorder, type = type)
  svdres <- svd(basiscormat, 1, 1)
  maxcor <- svdres$d[1]
  alphag <- as.vector(svdres$u)
  alphah <- as.vector(svdres$v)
  bexfun <- switch(type,
    legendre = udp::udplegendrebex,
    cosine = udp::udpcosinebex,
    stop("Unknown type of basis function")
  )
  Tg <- bexfun(alphag)
  # SVD determines alphag, alphah only up to a joint sign flip; orient them
  # so that g is increasing at u = 1. T = F(g) inherits g's sign of change
  # there because F is strictly increasing on the range of g, so
  # udp::udpderiv(Tg, 1) already carries the answer without needing g itself.
  if (udp::udpderiv(Tg, 1) < 0) {
    alphag <- -alphag
    alphah <- -alphah
    Tg <- bexfun(alphag)
  }
  methods::new("bex",
    type = type, maxcor = maxcor, alphag = alphag, alphah = alphah,
    Tg = Tg, Th = bexfun(alphah)
  )
}

#' Plot a Basis Expansion
#'
#' Draws the pair of expansion functions `g`, `h` over `[0, 1]` side by side,
#' or, with `udp = TRUE`, their \pkg{udp} transformations `Tg`, `Th`. Both
#' panels share a common y-axis range and are drawn in the default line
#' colour.
#'
#' @param x an object of class \linkS4class{bex}.
#' @param n number of grid points for plotting.
#' @param udp logical; plot the udp transformations `Tg`, `Th` rather than the raw expansions `g`, `h`.
#' @param xlab x-axis label, shared by both panels.
#' @param ylab y-axis label(s) for the `g`/`Tg` and `h`/`Th` panels
#'   respectively; a single value is recycled for both. Defaults to
#'   `c("g(u)", "h(u)")`, or `c("Tg(u)", "Th(u)")` when `udp = TRUE`.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' data <- copula::rCopula(1000, copula::tCopula(0, df = 1))
#' bex <- basisexpand(data, maxorder = 8, type = "legendre")
#' plot(bex)
#' plot(bex, udp = TRUE)
setMethod("plot", c(x = "bex", y = "missing"), function(x, n = 500L, udp = FALSE,
                                                         xlab = "u", ylab = NULL, ...) {
  if (is.null(ylab)) {
    ylab <- if (udp) c("Tg(u)", "Th(u)") else c("g(u)", "h(u)")
  }
  ylab <- rep_len(ylab, 2)
  u <- seq(0, 1, length.out = n)
  if (udp) {
    gvals <- udp::udptrans(x@Tg, u)
    hvals <- udp::udptrans(x@Th, u)
  } else {
    B <- sapply(seq_along(x@alphag), function(i) basisfunc(u, i, x@type))
    gvals <- as.vector(B %*% x@alphag)
    hvals <- as.vector(B %*% x@alphah)
  }
  # Tg, Th are udp transformations and so take values in exactly [0, 1];
  # range(gvals, hvals) would instead reflect their panel-spline
  # interpolation error (~1e-4), narrowing ylim fractionally inside [0, 1].
  ylim <- if (udp) c(0, 1) else range(gvals, hvals)
  # pty = "s" makes each side-by-side panel square in physical inches. Without
  # it, asp = 1 (below) forces 1 data-unit to mean the same physical distance
  # on both axes of an oblong panel, which it can only do by *stretching one
  # axis's displayed range past what ylim/xlim asked for -- visibly past 0 or
  # 1 for these udp transformations.
  op <- graphics::par(mfrow = c(1, 2), pty = "s")
  on.exit(graphics::par(op))
  if (udp) {
    plot(u, gvals, type = "l", xlab = xlab, ylab = ylab[1], ylim = ylim, asp = 1, xaxs = "i", yaxs = "i", main = "Tg", ...)
    plot(u, hvals, type = "l", xlab = xlab, ylab = ylab[2], ylim = ylim, asp = 1, xaxs = "i", yaxs = "i", main = "Th", ...)
  } else {
    plot(u, gvals, type = "l", xlab = xlab, ylab = ylab[1], ylim = ylim, xaxs = "i", yaxs = "i", main = "g", ...)
    plot(u, hvals, type = "l", xlab = xlab, ylab = ylab[2], ylim = ylim, xaxs = "i", yaxs = "i", main = "h", ...)
  }
})

#' @describeIn bex-class Apply the pair of udp transformations `Tg`, `Th` of
#'   a bex object to the two columns of a bivariate matrix, `cbind(udptrans(x@Tg,
#'   u[, 1]), udptrans(x@Th, u[, 2]))`. Unlike a single \link[udp]{udp-class}
#'   transformation, which acts on a vector, `u` here must be a two-column
#'   matrix since the two columns get different transformations. Use
#'   `udp::boundaryadjust()` on the result if values strictly inside `(0, 1)`
#'   are needed.
#' @param x an object of class \linkS4class{bex}.
#' @param u a two-column matrix with values in `(0, 1)`.
#' @importFrom udp udptrans
#' @export
setMethod("udptrans", "bex", function(x, u) {
  if (!is.matrix(u) || ncol(u) != 2) {
    stop("u must be a two-column matrix")
  }
  cbind(udp::udptrans(x@Tg, u[, 1]), udp::udptrans(x@Th, u[, 2]))
})
