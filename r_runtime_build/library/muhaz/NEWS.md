# muhaz 1.2.6.5

## Bug fixes

* Fixed a crash in `muhaz()` when a group has zero events. The default
  pilot bandwidth calculation (`bw.pilot <- endz/8/(nz^.2)`) divided by
  zero when `sum(delta) == 0`, producing `Inf` and crashing the
  downstream Fortran routine. `nz` is now clamped to a minimum of 1,
  which yields a finite bandwidth and a correctly flat zero hazard
  estimate when no events are present. Reported and fixed by
  Ken Beath <ken@kjbeath.id.au>.

## Documentation

* Fixed `man/muhaz.object.Rd`: replaced `\arguments` with `\value` /
  `\describe`, resolving an R CMD check NOTE ("Rd files without
  \usage... \arguments should not be documented without \usage").

## Package maintenance

* muhaz is now maintained by David Winsemius
  <dwinsemius@comcast.net>, following its orphaning by the previous
  maintainer.
* Added GitHub Actions / R-hub v2 CI configuration for cross-platform
  checks.

# muhaz 1.2.6.4

* Version on CRAN prior to this maintainer transition.
