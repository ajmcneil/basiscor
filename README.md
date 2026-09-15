# basiscor

<!-- badges: start -->
[![R-CMD-check](https://github.com/ajmcneil/basiscor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ajmcneil/basiscor/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`basiscor` computes **orthonormal basis correlations**: a family of
generalizations of Spearman's rank correlation, built from orthonormal basis
functions (shifted Legendre polynomials, or a cosine basis) of the two
uniform margins of a bivariate distribution. At degree 1 a basis correlation
is exactly Spearman's rho; higher degrees can pick up dependence -- including
non-monotonic dependence -- that ordinary correlation and Spearman's rho
miss entirely.

## Installation

Install the development version, and its dependency
[`udp`](https://github.com/ajmcneil/udp), from GitHub:

``` r
# install.packages("pak")
pak::pak("ajmcneil/basiscor")
```

## Example

`basiscor()` computes the population basis correlation between two degrees
`j`, `k` from a parametric copula object; `basiscordata()` computes its
sample analogue from data.

``` r
library(basiscor)

cop <- copula::claytonCopula(2)
basiscor(cop, 1, 1) # Spearman's rho
basiscor(cop, 3, 4) # a higher-degree basis correlation

X <- copula::rCopula(1000, cop)
basiscordata(X, 3, 4) # its sample analogue
```

`basiscormatrix()` computes a full matrix of basis correlations at once, and
`basisexpand()` finds the pair of basis-function expansions that maximizes
the generalized Spearman correlation between two margins -- useful for
uncovering non-monotonic dependence a sample's ordinary correlation misses
entirely. See `vignette("basis-expansions")` for a worked example, built on
the [`udp`](https://github.com/ajmcneil/udp) package's
uniform-distribution-preserving transformations.

## References

- McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
  of non-monotonic dependence. arXiv:2512.10828.
  <https://arxiv.org/abs/2512.10828>

## License

GPL-3 © Alexander J. McNeil, Johanna Nešlehová, Andrew D. Smith
