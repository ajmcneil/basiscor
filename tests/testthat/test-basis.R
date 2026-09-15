# sLegendre(), basisfunc(), basisderiv(), basisintegral()

test_that("sLegendre() matches known closed forms for low degree", {
  x <- seq(0, 1, length.out = 11)
  expect_equal(sLegendre(x, 0), rep(1, length(x)))
  expect_equal(sLegendre(x, 1), 2 * x - 1)
  expect_equal(sLegendre(x, 2), (3 * (2 * x - 1)^2 - 1) / 2)
})

test_that("sLegendre() derivative matches a central finite difference", {
  x <- seq(0.05, 0.95, length.out = 19)
  h <- 1e-6
  for (d in 1:6) {
    fd <- (sLegendre(x + h, d) - sLegendre(x - h, d)) / (2 * h)
    expect_equal(sLegendre(x, d, deriv = TRUE), fd, tolerance = 1e-4)
  }
})

test_that("sLegendre() rejects a negative degree", {
  expect_error(sLegendre(0.5, -1))
})

test_that("basisfunc() legendre matches sLegendre() scaled to orthonormal", {
  x <- seq(0, 1, length.out = 11)
  for (d in 0:6) {
    expect_equal(basisfunc(x, d, "legendre"), sLegendre(x, d) * sqrt(2 * d + 1))
  }
})

test_that("basisfunc() cosine matches its closed form", {
  x <- seq(0, 1, length.out = 11)
  for (d in 1:6) {
    expect_equal(
      basisfunc(x, d, "cosine"),
      cos(d * pi * x) * sqrt(2) * (-1)^d
    )
  }
})

test_that("basis functions are orthonormal on [0, 1]", {
  for (type in c("legendre", "cosine")) {
    for (jk in list(c(1, 1), c(2, 2), c(2, 3), c(1, 4))) {
      j <- jk[1]
      k <- jk[2]
      val <- stats::integrate(function(x) basisfunc(x, j, type) * basisfunc(x, k, type), 0, 1)$value
      expect_equal(val, as.numeric(j == k), tolerance = 1e-6)
    }
  }
})

test_that("basisderiv() matches a central finite difference of basisfunc()", {
  x <- seq(0.05, 0.95, length.out = 19)
  h <- 1e-6
  for (type in c("legendre", "cosine")) {
    for (d in 1:5) {
      fd <- (basisfunc(x + h, d, type) - basisfunc(x - h, d, type)) / (2 * h)
      expect_equal(basisderiv(x, d, type), fd, tolerance = 1e-4)
    }
  }
})

test_that("basisintegral() is the antiderivative of basisfunc() vanishing at 0", {
  x <- seq(0.05, 0.95, length.out = 19)
  h <- 1e-6
  for (type in c("legendre", "cosine")) {
    for (d in 1:5) {
      expect_equal(basisintegral(0, d, type), 0, tolerance = 1e-8)
      fd <- (basisintegral(x + h, d, type) - basisintegral(x - h, d, type)) / (2 * h)
      expect_equal(basisfunc(x, d, type), fd, tolerance = 1e-4)
    }
  }
})

test_that("basisintegral() over the full domain is 0 for degree >= 1 (orthogonality to P_0)", {
  for (type in c("legendre", "cosine")) {
    for (d in 1:5) {
      expect_equal(basisintegral(1, d, type), 0, tolerance = 1e-8)
    }
  }
})
