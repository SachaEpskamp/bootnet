# Regression test for audit item 7: nCores worker-process count.
#
# Fix in R/bootnet.R (and vendored R/parSim.R): nCores = k now starts k worker
# processes instead of k - 1. This exercises the parallel makeSOCKcluster path.
#
# NOTE: parallel bootstraps use independent RNG streams, so results are NOT
# expected to match a sequential run; we only assert the run completes and
# returns a valid "bootnet" object.

library(bootnet)

set.seed(1)
data <- matrix(rnorm(60 * 4), ncol = 4)

net <- estimateNetwork(data, default = "pcor", verbose = FALSE)
stopifnot(inherits(net, "bootnetResult"))

boot <- bootnet(
  net,
  nBoots = 4,
  nCores = 2,
  default = "pcor",
  verbose = FALSE
)

stopifnot(inherits(boot, "bootnet"))

cat("item07-ncores: OK\n")
