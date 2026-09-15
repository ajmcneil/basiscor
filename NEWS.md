# basiscor 1.1.5

* Rewritten on top of the **udp** package, which now supplies all of the
  uniform-distribution-preserving transformation machinery this package used
  to implement itself:
  * `distLegendre()`/`qLegendre()`/`legendre_initial()`/`legendre_roots()`
    and `basis.R`'s `udpfunc()`/`udpderiv()`/`udpstochinv()`/`cosine_roots()`
    are gone, replaced by `udp::udplegendre()`/`udpcosine()` and
    `udptrans()`/`udpinverse()`/`udpsi()`.
  * `rstochinv()`/`rstochinvbounds()` are gone -- now one-liners at the call
    site via `udp::udpsi()`.
  * `basisexpand()` now returns an S4 `bex` object (slots `type`, `maxcor`,
    `alphag`, `alphah`, `Tg`, `Th`) with `show()` and `plot()` methods, and a
    `udptrans()` method for applying the pair of transformations to a
    bivariate matrix at once. `udpbasis_sum()`, `basiscombo()`,
    `basiscombo_deriv()`, `bexfunc()`, `udpbex()`, `udpbexdata()` and
    `pbex()` are gone as a result.
  * `boundaryadjust()` has moved to the **udp** package (it isn't specific
    to basis correlations).
  * `extremalLegendre()` is back, now built on `udp::udpquantile()`.
  * The package no longer depends on **orthopolynom**; `sLegendre()` is a
    dependency-free three-term recurrence instead. It now depends on
    **udp** (`Remotes: ajmcneil/udp`, not yet on CRAN).
* `basiscor()`/`basiscorcopula()`/`basiscordata()` hardened:
  * `j`, `k` are validated as non-negative integers *before* coercion (a
    non-integer such as `2.5` used to be silently truncated), and the
    `j == 0`/`k == 0` orthonormality shortcut (implied by orthonormality:
    `1` if both are `0`, `0` if only one is) is now applied consistently to
    both the population and sample paths.
  * `type`/`method` are validated with `match.arg()` instead of silently
    falling through to Legendre, or erroring deep inside the integrand on
    an unmatched `method`.
  * The `tCopula` non-integer-degrees-of-freedom check (`pCopula` has no
    closed form there, so `basiscorcopula()` switches to `method = "d"`
    automatically) no longer assumes a fixed slot layout.
  * The `copula` parameter of `basiscorcopula()`/`basiscor_inner()`, which
    shadowed the **copula** package's own name, is renamed to `cop`.
* `basiscorcopula()`'s nested `stats::integrate()` is replaced by tensor
  Gauss-Legendre quadrature:
  * `method = "p"`: a fixed 40 nodes per axis, ~1.5-2x faster and accurate
    to within 6e-5 against the previous default-tolerance values (validated
    across several copula families).
  * `method = "d"`: nodes are chosen adaptively (starting at 20 per axis,
    doubling until successive estimates agree to `1e-4`, capped at 320),
    since how much resolution the density integrand needs depends on the
    copula's own dependence structure rather than predictably on `j`, `k`.
    The old code was not just slow here but unreliable -- it returned `NA`
    outright for some entirely ordinary cases. The new code never fails and
    is faster in every case tested where the old code succeeded (up to
    ~200x for the easiest cases, several-fold even for the hardest).
* `basiscordata()`'s default `method` is now `"T3"` (was `"T1"`): a
  correlation-type statistic (`stats::cor()`), rather than a raw mean
  product with no centring or scaling.
* `plot()` on a `bex` object now draws its two curves in separate,
  side-by-side panels in the default line colour, instead of superimposed
  in different colours -- both for the raw expansions and (`udp = TRUE`)
  the fitted udp transformations.
* Every plot with both axes on `[0, 1]` (in this package and in **udp**)
  now uses `asp = 1` so the unit square renders as an actual square.
* Documentation switched to Oxford (`-ize`) spelling, matching **udp**.
* Added a `basis-expansions` vignette illustrating the `bex` class, `plot()`
  with `udp = TRUE`, overlaying fitted against true udp transformations, and
  `udptrans()` on a fitted `bex` object, using the same folded-copula
  dataset as the **udp** package's own vignette.
* The `Calculations` vignette's Gumbel/Clayton discrepancy against
  `copula::rho()` is resolved: a large-sample Monte Carlo check shows
  `basiscor()`'s value is the more accurate of the two. It also gained
  worked examples of `basiscordata()`, `basiscormatrix()` and
  `extremalLegendre()`.
* General documentation improvements: fuller explanations for `basiscor()`,
  `basiscorcopula()`, `basiscordata()` and `basiscormatrix()`, a documented
  `method` argument for `basiscordata()` (`"T1"`-`"T6"`), and examples added
  to `basisfunc()`, `basisderiv()`, `basisintegral()` and `basiscordata()`
  where they were missing.

# basiscor 1.1.4

* Prior release.
