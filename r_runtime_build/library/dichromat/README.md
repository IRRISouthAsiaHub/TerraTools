# Simulating the Effects of Color Blindness in R

## Overview

The R package [dichromat](https://zeileis.codeberg.page/dichromat/) provides:

* Simulation of the effects of color blindness via function `dichromat()`
  for two types of red-green color blindness (`type = "deutan"` or `"protan"`)
  and green-blue color blindness (`type = "tritan"`).

* Data `dalton` from the work of Françoise Viénot and co-authors, underlying
  the `dichromat()` function. It provides a 256-color palette as it would
  appear with normal, protanope, deuteranope, and tritanope vision.

* A `colorschemes` list containing 17 color vectors suitable for people with
  deficient or anomalous red-green vision.

⚠️ **Important note:** The R package [colorspace](https://doi.org/10.18637/jss.v096.i01) provides better functions for simulating color vision deficiency: [simulate_cvd(), deutan(), protan(), tritan()](https://zeileis.codeberg.page/colorspace/reference/simulate_cvd.html). More modern robust color palettes are available in [hcl.colors() and palette.colors()](https://doi.org/10.32614/RJ-2023-071) in base R as well as further tools in colorspace.


## Installation

The stable version of `dichromat` is available on
[CRAN](https://CRAN.R-project.org/package=dichromat):

``` r
install.packages("dichromat")
```

The latest development version can be installed from
[R-universe](https://zeileis.R-universe.dev/dichromat):

``` r
install.packages("dichromat", repos = "https://zeileis.R-universe.dev")
```


## License

The package is available under the
[General Public License version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html).
