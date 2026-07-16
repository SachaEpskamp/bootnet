# CRAN comments for bootnet 1.9

This is a version update of the existing CRAN package bootnet. The current
CRAN release is 1.8. Version 1.9 is a major bug-fix release; Adela-Maria
Isvoranu has been added as a contributor in Authors@R.

## Release summary

* Fixed a crash on R >= 4.3 when the `missing` argument was supplied to
  `bootnet()` without an explicit `type` argument.
* Fixed `corStability()` for node-drop bootstraps (nNode/nNodes field
  mismatch produced NaN drop proportions).
* Fixed `plot(<case-drop bootnet>, plot = "interval")` (crash on a
  non-existent `nPeople` column; corrected x-axis limits).
* Fixed the internal `statTable()`: hybrid-failure fallback no longer
  overwrites valid rspbc rows; `bridgeInDegree`/`bridgeOutDegree` statistics
  are now actually stored.
* Fixed the parametric Ising bootstrap for data not encoded as 0/1: the
  response encoding is now detected from the data and passed to
  `IsingSampler()`; a new `responses` argument covers the manual parametric
  bootstrap. Results for 0/1-encoded data are unchanged.
* Fixed the graphicalVAR-datatype bootstrap paths (undefined `vars`
  variable; wrong grep-based subsetting of lagged design-matrix columns).
* `nCores` now uses the requested number of worker processes (previously one
  fewer was used).
* Cleaned up the vendored copy of parSim (renamed to `bootnet_parSim` so it
  no longer shadows the parSim package).
* Added a `maxErrors` argument to `bootnet()` and more informative
  bootstrap-failure errors/warnings.
* Added an extensive new test suite (14 test scripts), including
  backward-compatibility reference tests pinning bootnet 1.8 results for the
  EBICglasso, ggmModSelect, IsingFit, and mgm estimators as well as
  bootstrap outputs.

## Dependencies

The `IsingSampler (>= 0.2.3)` requirement is unchanged and still holds: the
parametric Ising bootstrap calls `IsingSampler()` with `n`, `graph`,
`thresholds` and `responses` only, all of which are present in the current
CRAN release of IsingSampler (0.2.4). The package has additionally been
checked against the pending IsingSampler 0.5.0 bug-fix release, whose
behaviour changes (the `method = "direct"` return type, and the honouring of
`constrain` and `beta`) do not affect any code path used by bootnet, which
calls `IsingSampler()` with the default `method = "MH"` and `EstimateIsing()`
with the default `beta = 1`.

## Test environments

* Local: macOS Sonoma 14.2.1 (aarch64-apple-darwin20), R 4.5.3 (2026-03-11)

## R CMD check --as-cran results

0 ERRORs, 0 WARNINGs, 1 NOTE.

The NOTE is a property of the local test machine rather than of the package:

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough
HTML Tidy.
```

macOS ships an outdated HTML Tidy, so this check is skipped locally; it is
expected to pass on the CRAN build machines.

All 14 test scripts and all examples pass.
