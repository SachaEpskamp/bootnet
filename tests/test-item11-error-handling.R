# Regression test for audit item 11: bootstrap error handling.
#
# Changes in R/bootnet.R (sequential and parallel bootstrap loops):
#   - New maxErrors argument (default 10, matching the previously hard-coded
#     retry limit per bootstrap sample).
#   - When the retry limit is reached, the stop() message now reports which
#     bootstrap failed and the underlying estimation error (previously an
#     uninformative "Maximum number of errors in bootstraps reached").
#   - One aggregate warning at the end of the run reporting the total number
#     of bootstrap estimations that failed and were resampled.
#
# The success path is unchanged (no RNG consumption added); that is covered
# by tests/test-backcompat.R.

library(bootnet)

set.seed(1)
data <- matrix(rnorm(50 * 4), ncol = 4)
colnames(data) <- paste0("V", 1:4)

# Custom estimator factory: fails deterministically on a given set of calls
# via a closure counter (sequential path only; closures do not carry state
# across parallel workers). Returns a simple correlation matrix as the
# weights matrix, which is all bootnet needs from a custom 'fun'.
makeFlaky <- function(failCalls) {
    count <- 0
    function(data) {
        count <<- count + 1
        if (count %in% failCalls) stop("flaky failure for testing")
        cor(data)
    }
}

## (a) Sequential: retries happen and succeed, and the aggregate warning
## reports the exact number of failures.
#
# Call accounting: bootnet(data, fun = ...) first estimates the sample
# network (call 1), so the failures must start at call 2 to land in the
# bootstrap phase. Calls 2:4 fail (3 retries within bootstrap 1), call 5
# succeeds, bootstraps 2 and 3 succeed (calls 6, 7).
set.seed(123)
warns <- character(0)
boot <- withCallingHandlers(
    bootnet(data, fun = makeFlaky(2:4), nBoots = 3, verbose = FALSE,
            maxErrors = 10),
    warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)
stopifnot(inherits(boot, "bootnet"))
stopifnot(length(boot$boots) == 3)
stopifnot(any(grepl("3 bootstrap estimation(s) failed and were resampled.",
                    warns, fixed = TRUE)))

## (b) Sequential: exceeding maxErrors stops with an informative error that
## includes the underlying estimation error message.
#
# Call 1 (sample estimation) succeeds; every bootstrap call fails. With
# maxErrors = 2 the first bootstrap gives up after 3 failed attempts.
set.seed(123)
err <- tryCatch(
    bootnet(data, fun = makeFlaky(2:1000), nBoots = 3, verbose = FALSE,
            maxErrors = 2),
    error = function(e) conditionMessage(e)
)
stopifnot(is.character(err))
stopifnot(grepl("Maximum number of retries", err, fixed = TRUE))
stopifnot(grepl("flaky failure for testing", err, fixed = TRUE))

## (c) Parallel: a healthy run still completes, bootResults keep their
## previous structure (list of bootnetResult objects, not wrapped lists),
## and no retry warning is emitted.
set.seed(42)
net <- estimateNetwork(data, default = "pcor", verbose = FALSE)
warnsPar <- character(0)
bootPar <- withCallingHandlers(
    bootnet(net, nBoots = 4, nCores = 2, default = "pcor", verbose = FALSE),
    warning = function(w) {
        warnsPar <<- c(warnsPar, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)
stopifnot(inherits(bootPar, "bootnet"))
stopifnot(length(bootPar$boots) == 4)
stopifnot(all(sapply(bootPar$boots, inherits, "bootnetResult")))
stopifnot(!any(grepl("failed and were resampled", warnsPar)))

cat("item11-error-handling: OK\n")
