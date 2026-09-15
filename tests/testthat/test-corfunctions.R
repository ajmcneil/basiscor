# validate_degrees(), basiscor(), basiscorcopula(), basiscordata(),
# basiscormatrix(), extremalLegendre(), discreteLegendre()

test_that("validate_degrees() rejects non-integer or negative j, k before coercion", {
  expect_error(validate_degrees(2.5, 3), "non-negative integers")
  expect_error(validate_degrees(-1, 0), "non-negative integers")
  expect_error(validate_degrees(2, -1), "non-negative integers")
  expect_error(validate_degrees(NA, 1), "non-negative integers")
})

test_that("validate_degrees() gives the orthonormality shortcut for j == 0 or k == 0", {
  expect_equal(validate_degrees(0, 0)$shortcut, 1)
  expect_equal(validate_degrees(0, 3)$shortcut, 0)
  expect_equal(validate_degrees(3, 0)$shortcut, 0)
  expect_true(is.na(validate_degrees(2, 3)$shortcut))
})

test_that("validate_degrees() coerces to integer", {
  dg <- validate_degrees(3, 4)
  expect_identical(dg$j, 3L)
  expect_identical(dg$k, 4L)
})

## basiscor() / basiscorcopula() -------------------------------------------

test_that("basiscor() dispatches on object type", {
  cop <- copula::claytonCopula(2)
  expect_equal(basiscor(cop, 2, 2), basiscorcopula(cop, 2, 2, copobj = TRUE))
  set.seed(1)
  X <- copula::rCopula(200, cop)
  expect_equal(basiscor(X, 2, 2), basiscordata(X, 2, 2))
  expect_error(basiscor(list(1, 2)), "Unknown object type")
})

test_that("basiscor() gives the orthonormality shortcut for degree 0", {
  cop <- copula::claytonCopula(2)
  expect_equal(basiscor(cop, 0, 0), 1)
  expect_equal(basiscor(cop, 0, 3), 0)
  expect_equal(basiscor(cop, 3, 0), 0)
})

test_that("basiscor() rejects non-integer/negative j, k even when k == 0", {
  cop <- copula::claytonCopula(2)
  expect_error(basiscor(cop, 2.5, 3))
  expect_error(basiscor(cop, -1, 0)) # used to silently return 0
})

test_that("basiscorcopula() validates method and type with match.arg()", {
  cop <- copula::claytonCopula(2)
  expect_error(basiscor(cop, 2, 3, type = "fourier"))
  expect_error(basiscor(cop, 2, 3, method = "x"))
  expect_error(basiscordata(matrix(runif(20), 10, 2), 2, 3, method = "T9"))
  expect_error(basiscordata(matrix(runif(20), 10, 2), 2, 3, type = "fourier"))
})

test_that("basiscor(cop, 1, 1) is Spearman's rho, matching a closed form", {
  rho <- 0.7
  cop <- copula::normalCopula(rho)
  expect_equal(basiscor(cop), (6 / pi) * asin(rho / 2), tolerance = 1e-4)
  expect_equal(basiscor(cop, method = "d"), (6 / pi) * asin(rho / 2), tolerance = 1e-4)
})

test_that("method = \"p\" and method = \"d\" agree closely", {
  for (cop in list(
    copula::claytonCopula(1.5), copula::gumbelCopula(1.5),
    copula::normalCopula(0.6), copula::tCopula(0.5, df = 4)
  )) {
    for (jk in list(c(1, 1), c(2, 3))) {
      vp <- basiscor(cop, jk[1], jk[2], method = "p")
      vd <- basiscor(cop, jk[1], jk[2], method = "d")
      expect_equal(vp, vd, tolerance = 5e-3)
    }
  }
})

test_that("basiscorcopula() switches a non-integer-df tCopula to method = \"d\" automatically", {
  cop <- copula::tCopula(0.5, df = 4.5)
  expect_equal(basiscor(cop, 2, 3), basiscor(cop, 2, 3, method = "d"))
})

test_that("basiscorcopula() requires a parCopula object when copobj = TRUE", {
  expect_error(
    basiscorcopula(matrix(1, 2, 2), 1, 1, copobj = TRUE),
    "parametric copula"
  )
})

test_that("basiscor() works with a user-supplied copula function (copobj = FALSE)", {
  # independence copula as a plain function
  indep <- function(u, v) u * v
  expect_equal(basiscor(indep, 1, 1), 0, tolerance = 1e-6)
  expect_equal(basiscor(indep, 2, 2), 0, tolerance = 1e-6)
})

test_that("basiscor() at high degree agrees between method p and d, cosine type too", {
  cop <- copula::frankCopula(3)
  expect_equal(
    basiscor(cop, 4, 5, type = "cosine", method = "p"),
    basiscor(cop, 4, 5, type = "cosine", method = "d"),
    tolerance = 1e-3
  )
})

## basiscordata() ------------------------------------------------------------

test_that("basiscordata() at j = k = 1, method = T3 equals sample Spearman correlation", {
  set.seed(1)
  X <- copula::rCopula(300, copula::claytonCopula(2))
  expect_equal(basiscordata(X, 1, 1), stats::cor(X[, 1], X[, 2], method = "spearman"))
})

test_that("basiscordata() requires a two-column matrix", {
  expect_error(basiscordata(data.frame(a = 1:5, b = 1:5), 1, 1), "matrix")
  expect_error(basiscordata(matrix(1:9, 3, 3), 1, 1), "2 columns")
})

test_that("basiscordata() gives the orthonormality shortcut for degree 0", {
  set.seed(1)
  X <- copula::rCopula(50, copula::claytonCopula(2))
  expect_equal(basiscordata(X, 0, 0), 1)
  expect_equal(basiscordata(X, 0, 3), 0)
})

test_that("basiscordata() T6 is legendre-only", {
  set.seed(1)
  X <- copula::rCopula(50, copula::claytonCopula(2))
  expect_error(basiscordata(X, 1, 1, type = "cosine", method = "T6"), "Legendre")
  expect_silent(basiscordata(X, 1, 1, type = "legendre", method = "T6"))
})

test_that("all six basiscordata() methods run and return a finite scalar", {
  set.seed(1)
  X <- copula::rCopula(100, copula::claytonCopula(2))
  for (m in c("T1", "T2", "T3", "T4", "T5", "T6")) {
    v <- basiscordata(X, 2, 3, method = m)
    expect_length(v, 1)
    expect_true(is.finite(v))
  }
})

## basiscormatrix() ----------------------------------------------------------

test_that("basiscormatrix() matches basiscor() entrywise", {
  cop <- copula::claytonCopula(2)
  m <- basiscormatrix(cop, maxorder = 3)
  for (j in 1:3) {
    for (k in 1:3) {
      expect_equal(m[j, k], basiscor(cop, j, k))
    }
  }
})

test_that("basiscormatrix(symmetric = TRUE) mirrors the upper triangle and halves the work", {
  cop <- copula::normalCopula(0.5) # exchangeable margins -> genuinely symmetric
  m <- basiscormatrix(cop, maxorder = 4, symmetric = TRUE)
  expect_equal(m, t(m))
  m_full <- basiscormatrix(cop, maxorder = 4, symmetric = FALSE)
  expect_equal(m, m_full, tolerance = 1e-6)
})

## extremalLegendre() ---------------------------------------------------------

test_that("extremalLegendre() gives max >= 0 >= min and max self-correlation is 1", {
  expect_equal(extremalLegendre(3, 3, "max"), 1, tolerance = 1e-3)
  bmax <- extremalLegendre(3, 4, "max")
  bmin <- extremalLegendre(3, 4, "min")
  expect_gt(bmax, 0)
  expect_lt(bmin, 0)
  expect_gt(bmax, bmin)
})

test_that("basiscor() values stay within the extremalLegendre() bounds", {
  bmax <- extremalLegendre(3, 4, "max")
  bmin <- extremalLegendre(3, 4, "min")
  for (cop in list(
    copula::claytonCopula(2), copula::gumbelCopula(3), copula::frankCopula(5),
    copula::normalCopula(0.7), copula::tCopula(0.6, df = 5)
  )) {
    v <- basiscor(cop, 3, 4)
    expect_true(v <= bmax + 1e-6 && v >= bmin - 1e-6)
  }
})

## discreteLegendre() / lff() -------------------------------------------------

test_that("discreteLegendre() matches a known reference value", {
  expect_equal(discreteLegendre(7, 20, 2), -21 / 57, tolerance = 1e-8)
})

test_that("lff() is the log of the falling factorial", {
  r <- 10
  k <- 3
  expect_equal(lff(r, k), log(r * (r - 1) * (r - 2)))
  expect_equal(lff(r, 0), 0)
})
