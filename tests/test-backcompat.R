# Backward-compatibility reference values captured from bootnet 1.8
# (commit 580f1fe) on R 4.6.1 — these results must not change unless a
# NEWS-documented, intentional break is made.
#
# All data are generated with fixed seeds; every reference value below was
# produced twice from the installed bootnet 1.8 and verified to be bitwise
# identical across runs (including mgm, which is seeded immediately before
# the call).

library(bootnet)

## ---------------------------------------------------------------------------
## Data generation (fixed seeds — do not change)
## ---------------------------------------------------------------------------

# Gaussian data: 6 variables, 200 rows, Toeplitz correlation matrix 0.6^|i-j|
genGauss <- function(){
  Sigma <- 0.6 ^ abs(outer(1:6, 1:6, "-"))
  set.seed(20260716)
  d <- as.data.frame(mvtnorm::rmvnorm(200, sigma = Sigma))
  colnames(d) <- paste0("V", 1:6)
  d
}

# Binary data: 5 nodes, chain Ising graph (weight 0.7), thresholds -0.7, n = 500
genBinary <- function(){
  ig <- matrix(0, 5, 5)
  ig[cbind(1:4, 2:5)] <- 0.7
  ig <- ig + t(ig)
  set.seed(31415)
  d <- as.data.frame(IsingSampler::IsingSampler(500, ig, rep(-0.7, 5)))
  colnames(d) <- paste0("I", 1:5)
  d
}

gaussData <- genGauss()
binData   <- genBinary()

## ---------------------------------------------------------------------------
## a. EBICglasso
## ---------------------------------------------------------------------------

net_ebic <- estimateNetwork(gaussData, default = "EBICglasso", verbose = FALSE)

ref_ebic_graph <- structure(c(0, 0.5013921896, 0, 0.0458909248, 0, 0, 0.5013921896,
0, 0.3989122116, 0.0509131063, 0.03146761, 0, 0, 0.3989122116,
0, 0.3440014768, 0, 0.0105654823, 0.0458909248, 0.0509131063,
0.3440014768, 0, 0.3829775915, 0.0133346708, 0, 0.03146761, 0,
0.3829775915, 0, 0.4493225882, 0, 0, 0.0105654823, 0.0133346708,
0.4493225882, 0), dim = c(6L, 6L))

stopifnot(isTRUE(all.equal(round(unname(net_ebic$graph), 10),
                           ref_ebic_graph, tolerance = 1e-8)))
cat("OK: EBICglasso graph\n")

## ---------------------------------------------------------------------------
## b. ggmModSelect
## ---------------------------------------------------------------------------

net_gms <- estimateNetwork(gaussData, default = "ggmModSelect", verbose = FALSE)

ref_gms_graph <- structure(c(0, 0.5588904338, 0, 0, 0, 0, 0.5588904338, 0, 0.4519451477,
0, 0, 0, 0, 0.4519451477, 0, 0.4105052623, 0, 0, 0, 0, 0.4105052623,
0, 0.4294993889, 0, 0, 0, 0, 0.4294993889, 0, 0.5029608362, 0,
0, 0, 0, 0.5029608362, 0), dim = c(6L, 6L))

stopifnot(isTRUE(all.equal(round(unname(net_gms$graph), 10),
                           ref_gms_graph, tolerance = 1e-8)))
cat("OK: ggmModSelect graph\n")

## ---------------------------------------------------------------------------
## c. IsingFit
## ---------------------------------------------------------------------------

net_ising <- suppressWarnings(
  estimateNetwork(binData, default = "IsingFit", verbose = FALSE))

ref_ising_graph <- structure(c(0, 0.8472137518, 0, 0, 0, 0.8472137518, 0, 0.3789165163,
0, 0, 0, 0.3789165163, 0, 0.6413969983, 0, 0, 0, 0.6413969983,
0, 0.4259031921, 0, 0, 0, 0.4259031921, 0), dim = c(5L, 5L))

ref_ising_intercepts <- c(-0.7827314985, -0.7543547939, -0.4739905321,
                          -0.7296237443, -0.5371955831)

stopifnot(isTRUE(all.equal(round(unname(net_ising$graph), 10),
                           ref_ising_graph, tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(round(unname(unlist(net_ising$intercepts)), 10),
                           ref_ising_intercepts, tolerance = 1e-8)))
cat("OK: IsingFit graph + intercepts\n")

## ---------------------------------------------------------------------------
## d. mgm (seed set immediately before the call; criterion = EBIC default)
## ---------------------------------------------------------------------------

set.seed(42)
net_mgm <- suppressWarnings(suppressMessages(
  estimateNetwork(binData, default = "mgm", verbose = FALSE)))

ref_mgm_graph <- structure(c(0, 0.423607139, 0, 0, 0, 0.423607139, 0, 0.1894588907,
0, 0, 0, 0.1894588907, 0, 0.3206992763, 0, 0, 0, 0.3206992763,
0, 0.2129518263, 0, 0, 0, 0.2129518263, 0), dim = c(5L, 5L))

stopifnot(isTRUE(all.equal(round(unname(net_mgm$graph), 10),
                           ref_mgm_graph, tolerance = 1e-8)))
cat("OK: mgm graph\n")

## ---------------------------------------------------------------------------
## e. Bootstrap results (nonparametric + person/case-drop) on EBICglasso net
## ---------------------------------------------------------------------------

set.seed(1234)
boot_np <- suppressMessages(
  bootnet(net_ebic, nBoots = 6, type = "nonparametric",
          statistics = c("edge", "strength"), verbose = FALSE))

# Sorted bootTable values (sorting makes this robust to row ordering)
ref_np_values <- c(-0.13939805, -0.10966104, -0.06282389, -0.03739441, -0.03564222,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00301201, 0.0088298, 0.00883062,
0.01299946, 0.01586979, 0.01590527, 0.01954646, 0.03346298, 0.03604092,
0.03668636, 0.03819781, 0.03957764, 0.03986875, 0.04030457, 0.04731852,
0.05130405, 0.0623849, 0.06917367, 0.09438307, 0.10960173, 0.1112265,
0.11175038, 0.11859931, 0.11885293, 0.1198825, 0.24607124, 0.28089999,
0.30633193, 0.33044699, 0.34043903, 0.34078203, 0.35364731, 0.37190248,
0.37300506, 0.37776775, 0.38969664, 0.40522278, 0.4107018, 0.41833464,
0.41899725, 0.42136412, 0.42919942, 0.43874591, 0.46995338, 0.48533852,
0.48840078, 0.4921327, 0.49328085, 0.49478216, 0.49496943, 0.51431166,
0.51918522, 0.52293325, 0.52983258, 0.54763095, 0.54901693, 0.56858632,
0.58465914, 0.5883564, 0.59526182, 0.59548093, 0.59597349, 0.6803367,
0.69123116, 0.69573331, 0.69853425, 0.72050956, 0.75085337, 0.75088232,
0.76085628, 0.76180316, 0.77848275, 0.79152637, 0.86249941, 0.8680938,
0.86872975, 0.87372246, 0.87974868, 0.90398624, 0.90438782, 0.92202198,
0.93543972, 0.93566801, 0.94559053, 0.97514951, 0.9868398, 1.0023142,
1.01596065, 1.06289669, 1.08311673, 1.14571073)

ref_np_dim  <- c(126L, 12L)
ref_np_cols <- c("name", "type", "node1", "node2", "value", "id", "nNode",
                 "nPerson", "rank_avg", "rank_min", "rank_max", "graph")

stopifnot(isTRUE(all.equal(round(sort(boot_np$bootTable$value), 8),
                           ref_np_values, tolerance = 1e-8)))
stopifnot(identical(dim(boot_np$bootTable), ref_np_dim))
stopifnot(identical(colnames(boot_np$bootTable), ref_np_cols))
cat("OK: nonparametric bootstrap values + bootTable structure\n")

set.seed(5678)
boot_cd <- suppressMessages(
  bootnet(net_ebic, nBoots = 6, type = "person", caseN = 3,
          verbose = FALSE))

ref_cs <- c(edge = 0.75, strength = 0.75)

cs <- round(suppressWarnings(corStability(boot_cd, verbose = FALSE)), 8)
stopifnot(isTRUE(all.equal(cs[names(ref_cs)], ref_cs, tolerance = 1e-8)))
cat("OK: person-drop bootstrap CS-coefficients\n")

cat("All backward-compatibility reference checks passed.\n")
