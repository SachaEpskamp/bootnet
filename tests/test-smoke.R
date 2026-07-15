# Fast smoke test: bootnet runs end-to-end with a trivial pcor network,
# summary() returns a data frame, and plot() returns a buildable ggplot.

library(bootnet)

# Small seeded Gaussian dataset: 5 variables, 100 rows
Sigma <- 0.5 ^ abs(outer(1:5, 1:5, "-"))
set.seed(20260716)
smokeData <- as.data.frame(mvtnorm::rmvnorm(100, sigma = Sigma))
colnames(smokeData) <- paste0("V", 1:5)

# Estimate:
net <- estimateNetwork(smokeData, default = "pcor", verbose = FALSE)
stopifnot(is.matrix(net$graph), identical(dim(net$graph), c(5L, 5L)))
cat("OK: estimateNetwork (pcor) runs\n")

# Bootstrap:
set.seed(1)
boots <- suppressMessages(
  bootnet(net, nBoots = 3, type = "nonparametric", verbose = FALSE))
stopifnot(inherits(boots, "bootnet"),
          is.data.frame(as.data.frame(boots$bootTable)),
          nrow(boots$bootTable) > 0)
cat("OK: bootnet (nonparametric, nBoots = 3) runs\n")

# summary() returns a data frame:
smry <- summary(boots)
stopifnot(is.data.frame(as.data.frame(smry)), nrow(smry) > 0)
cat("OK: summary() returns a data frame\n")

# plot() returns a ggplot object that builds without error:
plt <- plot(boots, order = "sample")
stopifnot(inherits(plt, "ggplot"))
built <- ggplot2::ggplot_build(plt)
stopifnot(!is.null(built$data), length(built$data) > 0)
cat("OK: plot() returns a buildable ggplot\n")

cat("Smoke test passed.\n")
