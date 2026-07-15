# Regression test for audit item 9: internal cleanups (behavior-neutral).
#
# Covers:
#  - R/zzz.R no longer creates NULL bindings shadowing base/stats functions
#    (sd, quantile, runif, na.omit, data) in the bootnet namespace.
#  - estimateNetwork(default = "pcor") still runs (backward-compat detail;
#    exact values are asserted in test-backcompat.R).
#  - bootnet(<raw data>, default = "graphicalVAR") errors with the informative
#    message, not match.arg()'s generic "should be one of" error.

library(bootnet)

## --- No NULL shadow bindings remain in the namespace --------------------------
ns <- asNamespace("bootnet")
for (nm in c("sd", "quantile", "runif", "na.omit", "data")) {
  # Either the name is not bound in the package namespace at all (it resolves
  # to the imported stats/utils function via the imports environment), or if it
  # is bound it must be the actual function -- never a NULL shadow.
  ok <- !exists(nm, envir = ns, inherits = FALSE) ||
    is.function(get(nm, envir = ns, inherits = FALSE))
  stopifnot(ok)
}
# The old fooling block explicitly assigned NULL; make sure none of these are NULL.
for (nm in c("sd", "quantile", "runif", "na.omit", "data", "P", "fill", "label")) {
  if (exists(nm, envir = ns, inherits = FALSE)) {
    stopifnot(!is.null(get(nm, envir = ns, inherits = FALSE)))
  }
}

## --- pcor estimation smoke check ---------------------------------------------
set.seed(1)
d <- matrix(rnorm(60 * 4), ncol = 4)
net <- estimateNetwork(d, default = "pcor", verbose = FALSE)
stopifnot(inherits(net, "bootnetResult"))
stopifnot(is.matrix(net$graph))

## --- graphicalVAR informative error before match.arg -------------------------
err <- tryCatch(
  bootnet(d, default = "graphicalVAR"),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err))
stopifnot(grepl("only supported for output of estimateNetwork", err, fixed = TRUE))
# Confirm it is NOT the generic match.arg error.
stopifnot(!grepl("should be one of", err))

cat("item09-cleanups: OK\n")
