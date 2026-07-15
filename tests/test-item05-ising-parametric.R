# Test for audit item 5: parametric Ising bootstrap must respect the response
# encoding of the data (e.g., -1/+1), the 0/1 path must remain bit-identical,
# and the parametric GGM bootstrap must reject correlation networks.
library(bootnet)

## ---------------------------------------------------------------------------
## 1. -1/+1 encoded Ising data: parametric bootstrap samples must be drawn on
##    the same encoding as the sample network parameterization.
##
## Before the fix, IsingSampler() was always called with its default
## responses = c(0L, 1L), simulating 0/1 data from the -1/+1 parameterized
## graph/thresholds: a much weaker model on that scale, so eLasso found almost
## no edges in the bootstrap samples (1/18 nonzero in this design).
## ---------------------------------------------------------------------------
igPM <- matrix(0, 4, 4)
for (i in 1:3) igPM[i, i + 1] <- igPM[i + 1, i] <- 0.4
set.seed(42)
pmData <- as.data.frame(
  IsingSampler::IsingSampler(500, igPM, rep(0, 4), responses = c(-1L, 1L)))
names(pmData) <- paste0("I", 1:4)
stopifnot(identical(sort(unique(unlist(pmData))), c(-1L, 1L)))

netPM <- suppressWarnings(
  estimateNetwork(pmData, default = "IsingFit", verbose = FALSE))
sampEdges <- netPM$graph[upper.tri(netPM$graph)]
stopifnot(any(sampEdges != 0))

set.seed(43)
btPM <- suppressWarnings(suppressMessages(
  bootnet(netPM, type = "parametric", nBoots = 3, verbose = FALSE)))

bootEdges <- btPM$bootTable$value[btPM$bootTable$type == "edge"]

# The bootstrap networks must recover edges on the sample network's scale:
# at least two thirds of the true (chain) edge slots nonzero across boots...
stopifnot(sum(bootEdges != 0) >= 6)

# ...and the median absolute nonzero edge weight on the same scale as the
# sample network (ratio within [0.25, 4]):
sampScale <- median(abs(sampEdges[sampEdges != 0]))
bootScale <- median(abs(bootEdges[bootEdges != 0]))
ratio <- bootScale / sampScale
stopifnot(is.finite(ratio), ratio >= 0.25, ratio <= 4)
cat("OK: -1/+1 parametric Ising bootstrap on sample-network scale (ratio =",
    round(ratio, 3), ")\n")

## ---------------------------------------------------------------------------
## 2. 0/1 regression guard: the parametric bootstrap for 0/1 data must be
##    bit-identical to before the fix. Reference values captured on commit
##    2a2680f (pre-change), R 4.6.0, IsingSampler 0.5.0.
## ---------------------------------------------------------------------------
ig01 <- matrix(0, 5, 5)
for (i in 1:4) ig01[i, i + 1] <- ig01[i + 1, i] <- 0.7
set.seed(2025)
binData <- as.data.frame(IsingSampler::IsingSampler(500, ig01, rep(-0.7, 5)))
names(binData) <- paste0("I", 1:5)

net01 <- suppressWarnings(
  estimateNetwork(binData, default = "IsingFit", verbose = FALSE))

set.seed(123)
bt01 <- suppressWarnings(suppressMessages(
  bootnet(net01, type = "parametric", nBoots = 3, verbose = FALSE)))

tab01 <- bt01$bootTable
tab01 <- tab01[order(tab01$name, tab01$type, tab01$node1, tab01$node2), ]

ref01 <- c(0, 0, 0, 0, 0.5463086932, 0, 0, 0.7820502957, 0, 0, -0.7263327295,
           -0.4401839977, -0.7852529674, -0.7762135698, -0.41380544, 0,
           0.5463086932, 1.3283589889, 0.7820502957, 0, 0, 0, 0, 0, 0.7637351289,
           0, 0, 0.6152623018, 0, 0, -0.506561225, -0.5482472102, -0.9065888822,
           -0.4969410404, -0.447312218, 0, 0.7637351289, 1.3789974307, 0.6152623018,
           0, 0, 0, 0, 0, 1.0813972727, 0, 0, 0.4452806869, 0, 0, -0.5753641449,
           -0.6736096209, -0.9968809715, -0.4265640625, -0.2411620568, 0,
           1.0813972727, 1.5266779596, 0.4452806869, 0)

stopifnot(nrow(tab01) == length(ref01))
stopifnot(isTRUE(all.equal(round(tab01$value, 10), ref01, tolerance = 1e-8)))
cat("OK: 0/1 parametric Ising bootstrap identical to pre-change reference\n")

## ---------------------------------------------------------------------------
## 3. GGM guard: parametric bootstrap must reject correlation networks.
## ---------------------------------------------------------------------------
set.seed(7)
gausData <- as.data.frame(
  mvtnorm::rmvnorm(200, sigma = diag(0.5, 4) + matrix(0.5, 4, 4)))
names(gausData) <- paste0("V", 1:4)

netCor <- estimateNetwork(gausData, default = "cor", verbose = FALSE)

err <- tryCatch({
  suppressWarnings(suppressMessages(
    bootnet(netCor, type = "parametric", nBoots = 2, verbose = FALSE)))
  NULL
}, error = function(e) conditionMessage(e))

stopifnot(!is.null(err))
stopifnot(grepl("requires a partial-correlation network", err, fixed = TRUE))
stopifnot(grepl("default = 'cor'", err, fixed = TRUE))
cat("OK: parametric GGM bootstrap rejects correlation networks\n")

cat("All item 5 checks passed.\n")
