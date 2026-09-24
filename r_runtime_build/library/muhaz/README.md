# muhaz

Smooth hazard function estimation from right-censored survival data.

## Overview

`muhaz` estimates the hazard function from right-censored data using
kernel-based methods, implementing the bandwidth selection algorithms and
boundary kernel formulations described in Mueller and Wang (1994). Options
include:

- Three bandwidth methods: local MSE-optimal, global IMSE-optimal, and
  k-nearest-neighbor
- Three boundary correction types: none, left only, or both ends
- Four kernel shapes: rectangle, Epanechnikov, biquadratic, triquadratic

A complementary set of piecewise-exponential estimators (`pehaz`,
`plot.pehaz`) is also provided for quick exploratory comparison.

Original S code by Kenneth Hess (M.D. Anderson Cancer Center); R port by
R. Gentleman. Currently maintained by David Winsemius.

## Installation

Install the released version from CRAN:

```r
install.packages("muhaz")
```

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("dwinsemius/muhaz")
```

## Quick start

```r
library(muhaz)
data(cancer, package = "survival")

# Locally optimal bandwidth (default)
fit <- muhaz(ovarian$futime, ovarian$fustat)
plot(fit)
summary(fit)

# Globally optimal bandwidth
fit_global <- muhaz(ovarian$futime, ovarian$fustat, bw.method = "g")

# Fixed bandwidth
fit_fixed <- muhaz(ovarian$futime, ovarian$fustat, bw.method = "g", bw.grid = 5)
```

## References

1. Mueller HG, Wang JL. Hazard rates estimation under random censoring with
   varying kernels and bandwidths. *Biometrics* 1994; **50**: 61--76.

2. Gefeller O, Dette H. Nearest neighbour kernel estimation of the hazard
   function from censored data. *J Statist Comput Simul* 1992; **43**: 93--101.

3. Hess KR, Serachitopol DM, Brown BW. Hazard function estimators: a
   simulation study. *Statistics in Medicine* 1999.

## License

GPL. See `COPYING` for details.
