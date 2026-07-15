# Regression test for audit item 1:
# bootnet() crashed on R >= 4.3 when `missing` was supplied without an explicit
# `type` argument, because the missing-handling check ran before match.arg(type),
# so `type != "parametric"` operated on the full choices vector and raised
# "length = 6 in coercion to logical(1)".

library(bootnet)

set.seed(1)
data <- as.data.frame(matrix(rnorm(50 * 4), nrow = 50, ncol = 4))
colnames(data) <- paste0("V", 1:4)

# Call bootnet() with `missing` supplied but NO `type` argument.
res <- tryCatch(
  bootnet(data, default = "pcor", missing = "listwise", nBoots = 2, verbose = FALSE),
  error = function(e) e
)

# Must not raise the coercion error, and must actually succeed.
if (inherits(res, "error")) {
  stop("bootnet() raised an error: ", conditionMessage(res))
}

stopifnot(inherits(res, "bootnet"))

cat("test-item01-missing-arg.R passed\n")
