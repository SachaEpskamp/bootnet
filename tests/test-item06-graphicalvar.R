# Regression tests for audit item 6: graphicalVAR-datatype bootstrap paths.
#
# Bugs fixed in R/bootnet.R:
#   1. Node-drop setup used `length(vars)` with `vars` undefined in bootnet()'s
#      scope -> "object 'vars' not found" for any node-drop bootstrap on
#      graphicalVAR data.
#   2. Node-drop subsetting of data_l used grep() with a character VECTOR as
#      pattern, so only the first sampled node's lag columns were kept ->
#      wrong design matrix.
#   3. The tsData 'vars' field was not subset, but graphicalVAR::graphicalVAR()
#      uses it to locate the lag-1 beta columns -> subscript error.

library(bootnet)

if (!requireNamespace("graphicalVAR", quietly = TRUE)) {
  message("Package 'graphicalVAR' not available; skipping item 6 tests.")
} else {

  set.seed(42)
  Mod <- graphicalVAR::randomGVARmodel(4, probKappaEdge = 0.3, probBetaEdge = 0.3)
  simData <- graphicalVAR::graphicalVARsim(100, Mod$beta, Mod$kappa)

  net <- suppressWarnings(estimateNetwork(
    simData,
    default = "graphicalVAR",
    nLambda = 4,
    verbose = FALSE
  ))
  stopifnot(inherits(net, "bootnetResult"))
  stopifnot(net$datatype == "graphicalVAR")

  # Sanity check of the stored tsData-like object and its lag-column naming
  # (intercept column "1" plus "<var>_lag1" per variable):
  stopifnot(!is.null(net$data$vars))
  stopifnot(identical(
    names(net$data$data_l),
    c("1", paste0(net$data$vars, "_lag1"))
  ))

  ## 1) Nonparametric bootstrap completes:
  set.seed(1)
  bNP <- suppressWarnings(bootnet(
    net, nBoots = 2, type = "nonparametric",
    verbose = FALSE, statistics = "edge"
  ))
  stopifnot(inherits(bNP, "bootnet"))
  stopifnot(length(bNP$boots) == 2)

  ## 2) Node-drop bootstrap completes (previously: "object 'vars' not found"):
  set.seed(2)
  bND <- suppressWarnings(bootnet(
    net, nBoots = 2, type = "node",
    verbose = FALSE, statistics = "edge",
    memorysaver = FALSE  # keep boot data so retained columns can be inspected
  ))
  stopifnot(inherits(bND, "bootnet"))
  stopifnot(length(bND$boots) == 2)

  # Each node-drop bootstrap must have retained exactly the sampled nodes'
  # columns: contemporaneous data, intercept + lag columns, vars field, and
  # graph dimensions all matching the sampled node set:
  for (i in seq_along(bND$boots)) {
    boot_i <- bND$boots[[i]]
    sampled <- boot_i$labels
    stopifnot(all(sampled %in% net$labels))
    stopifnot(boot_i$nNode == length(sampled))
    stopifnot(boot_i$nNode >= 2, boot_i$nNode <= 3) # subNodes = 2:(p-1)

    # Graph dimensions match the sampled node count:
    stopifnot(identical(dim(boot_i$graph$contemporaneous),
                        c(length(sampled), length(sampled))))
    stopifnot(identical(dim(boot_i$graph$temporal),
                        c(length(sampled), length(sampled))))

    # The bootstrapped tsData kept exactly the sampled nodes' columns:
    bd <- boot_i$data
    stopifnot(identical(bd$vars, sampled))
    stopifnot(identical(colnames(bd$data_c), sampled))
    stopifnot(identical(names(bd$data_l), c("1", paste0(sampled, "_lag1"))))
  }

  message("test-item06-graphicalvar.R: all tests passed.")
}
