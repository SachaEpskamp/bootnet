# bootnet 1.9

This is a version update of the existing CRAN package bootnet, maintained by
Sacha Epskamp <mail@sachaepskamp.com>. Version 1.9 is a major bug-fix release;
Adela Maria Isvoranu has been added as a contributor in Authors@R.

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

The `IsingSampler (>= 0.2.3)` requirement is unchanged. The package was
tested against CRAN IsingSampler 0.4.0's API and additionally against the
pending IsingSampler 0.5.0 bug-fix release.

## Test environments

* Local: macOS 26.4.1 (Apple Silicon), R 4.6.1 (2026-06-24)
* Checked with `R CMD check --as-cran --no-manual`
  (`--no-manual` because pdflatex is not installed locally, so the PDF
  manual could not be built; the Rd files pass checkRd without warnings)

## R CMD check results

Status: 1 NOTE (0 ERRORs, 0 WARNINGs)

* checking top-level files ... NOTE
  "Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed."

  This is a local-toolchain NOTE only: pandoc is not installed on the check
  machine, so README.md could not be validated locally. It does not indicate
  a problem with the package and will not occur on CRAN's builders, which
  have pandoc available.
