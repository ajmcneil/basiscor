# bex class, basisexpand(), plot(bex), udptrans(bex, .), show(bex)

test_that("basisexpand() returns a bex object with the expected slots", {
  set.seed(1)
  X <- copula::rCopula(500, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 5)
  expect_s4_class(bex, "bex")
  expect_identical(bex@type, "legendre")
  expect_length(bex@alphag, 5)
  expect_length(bex@alphah, 5)
  expect_s4_class(bex@Tg, "udplegendrebex")
  expect_s4_class(bex@Th, "udplegendrebex")
  expect_true(bex@maxcor >= 0 && bex@maxcor <= 1 + 1e-8)
})

test_that("basisexpand() respects type = \"cosine\"", {
  set.seed(1)
  X <- copula::rCopula(300, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 4, type = "cosine")
  expect_identical(bex@type, "cosine")
  expect_s4_class(bex@Tg, "udpcosinebex")
  expect_s4_class(bex@Th, "udpcosinebex")
})

test_that("basisexpand()'s maxcor is at least the plain Spearman correlation", {
  set.seed(1)
  X <- copula::rCopula(500, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 6)
  rho_s <- abs(stats::cor(X[, 1], X[, 2], method = "spearman"))
  expect_gte(bex@maxcor, rho_s - 1e-8)
})

test_that("basisexpand()'s sign convention makes g increasing at u = 1", {
  set.seed(1)
  for (cop in list(copula::claytonCopula(2), copula::gumbelCopula(3), copula::frankCopula(-3))) {
    X <- copula::rCopula(400, cop)
    bex <- basisexpand(X, maxorder = 6)
    expect_gte(udp::udpderiv(bex@Tg, 1), 0)
  }
})

test_that("basisexpand() recovers a known non-monotonic folding transform", {
  skip_on_cran()
  T5 <- udp::udplegendre(5)
  T3 <- udp::udpcosine(6)
  set.seed(1)
  W <- copula::rCopula(1000, copula::gumbelCopula(3))
  X <- cbind(udp::udpsi(T5, W[, 1]), udp::udpsi(T3, W[, 2]))
  bex <- basisexpand(X) # default maxorder = 8

  u <- seq(0, 1, length.out = 400)
  cor_match <- stats::cor(udp::udptrans(bex@Tg, u), udp::udptrans(T5, u))
  cor_cross <- stats::cor(udp::udptrans(bex@Tg, u), udp::udptrans(T3, u))
  expect_gt(abs(cor_match), 0.9)
  expect_lt(abs(cor_cross), 0.2)

  # the fitted transform pair should approach the copula's own correlation
  expect_gt(bex@maxcor, 0.7)
})

## udptrans(bex, .) -----------------------------------------------------------

test_that("udptrans() on a bex object applies Tg, Th to the two columns", {
  set.seed(1)
  X <- copula::rCopula(200, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 4)
  V <- udp::udptrans(bex, X)
  expect_equal(dim(V), dim(X))
  expect_equal(V[, 1], udp::udptrans(bex@Tg, X[, 1]))
  expect_equal(V[, 2], udp::udptrans(bex@Th, X[, 2]))
  expect_true(all(V >= 0 & V <= 1))
})

test_that("udptrans() on a bex object requires a two-column matrix", {
  set.seed(1)
  X <- copula::rCopula(200, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 4)
  expect_error(udp::udptrans(bex, X[, 1]), "two-column matrix")
  expect_error(udp::udptrans(bex, matrix(1, 3, 3)), "two-column matrix")
})

## show() / plot() -------------------------------------------------------------

test_that("show(bex) prints without error and mentions the key fields", {
  set.seed(1)
  X <- copula::rCopula(200, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 4)
  out <- capture.output(print(bex))
  expect_true(any(grepl("bex", out)))
  expect_true(any(grepl("maxcor", out)))
})

test_that("plot(bex) runs for both the raw and udp = TRUE panels", {
  set.seed(1)
  X <- copula::rCopula(200, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 4)
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  })
  expect_no_error(plot(bex))
  expect_no_error(plot(bex, udp = TRUE))
})

test_that("plot(bex, udp = TRUE) keeps both axes within [0, 1]", {
  # Tg, Th take values in exactly [0, 1]; a side-by-side (non-square) panel
  # under asp = 1 can stretch one axis's displayed range past the requested
  # ylim/xlim if the panel isn't forced square -- verify par("usr") (the
  # actually rendered extent of the last panel drawn, Th) stays inside
  # [0, 1] up to floating-point noise. Both panels share the same ylim/xlim
  # and pty = "s" setting, so this covers Tg's panel equally.
  set.seed(1)
  X <- copula::rCopula(200, copula::claytonCopula(2))
  bex <- basisexpand(X, maxorder = 4)
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  })
  plot(bex, udp = TRUE)
  usr <- graphics::par("usr")
  expect_true(all(usr > -1e-6 & usr < 1 + 1e-6))
})
