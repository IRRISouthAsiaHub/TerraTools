# dichromat 2.0-1

* Added pointers to the `colorspace` package and base R in a new
  startup message. The `colorspace` package provides better functions
  for simulating color vision deficiency: `simulate_cvd()`, `deutan()`,
  `protan()`, `tritan()`. More modern robust color palettes are available
  in `hcl.colors()` and `palette.colors()` in base R as well as further
  tools in `colorspace`. References providing more details:
  <https://doi.org/10.18637/jss.v096.i01> and
  <https://doi.org/10.32614/RJ-2023-071>.

* Code repository changed from R-Forge to Codeberg at:
  <https://codeberg.org/zeileis/dichromat/>

* Added basic altdoc page with overview and documentation at:
  <https://zeileis.codeberg.page/dichromat/>

* Achim Zeileis takes over maintenance.
