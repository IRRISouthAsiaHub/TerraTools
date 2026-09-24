<!-- README.md is generated from README.Rmd. Please edit that file -->

# `robust`: Port of the S+ “Robust Library”

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/robust)](https://cran.r-project.org/package=robust)
[![R-CMD-check](https://github.com/valentint/robust/workflows/R-CMD-check/badge.svg)](https://github.com/valentint/robust/actions)
[![downloads](https://cranlogs.r-pkg.org/badges/robust)](https://cran.r-project.org/package=robust)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/robust)](https://cran.r-project.org/package=robust)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

This package contains the Robust Library version 0.4

-   Contributors:
    -   Jeff Wang <jwang@statsci.com>
    -   Ruben Zamar <ruben@stat.ubc.ca>
    -   Alfio Marazzi <Alfio.Marazzi@inst.hospvd.ch>
    -   Victor Yohai <vyohai@dm.uba.ar>
    -   Matias Salibian-Barrera <matias@stat.ubc.ca>
    -   Ricardo Maronna <maron@mate.unlp.edu.ar>
    -   Eric Zivot <ezivot@u.washington.edu>
    -   David Rocke <dmrocke@ucdavis.edu>
    -   Doug Martin <doug@statsci.com>
    -   Kjell Konis <kjell.konis@icloud.com>

------------------------------------------------------------------------

-   This package contains the following robust methods:
    -   Robust Covariance estimation (scatter and location)
    -   Robust Linear Regression
    -   Robust Generalized Linear Models
    -   Robust Gamma, Weibull, and lognormal parameter estimation
-   A method for side-by-side comparison of robust and classical models
    is also provided. Please see Robust.pdf for further details.

## Installation

The `robust` package is on CRAN (The Comprehensive R Archive Network)
and the latest release can be easily installed using the command

    install.packages("robust")

## Building from source

To install the latest stable development version from GitHub, you can
pull this repository and install it using

    ## install.packages("remotes")
    remotes::install_github("valentint/robust")

Of course, if you have already installed `remotes`, you can skip the
first line (I have commented it out).

## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for
additional features, please submit an issue via the
[*Issues*](https://github.com/valentint/robust/issues) tab of this
repository. Please have a look at existing issues first to see if your
problem or feature request has already been discussed.

### Contribute to the package

If you want to contribute to the package, you can fork this repository
and create a pull request after implementing the desired functionality.

### Ask for help

If you need help using the package, or if you are interested in
collaborations related to this project, please get in touch with the
package maintainer.
