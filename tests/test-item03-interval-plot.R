# Audit item 3: plot(<case-drop bootnet>, plot = "interval") crashed with
# "Column `nPeople` doesn't exist" because R/plotMethod.R referenced a column
# never produced by summary.bootnet / statTable (actual name is nPerson).
# The x-axis break limits were also based on the number of nodes rather than
# the number of persons. See R/plotMethod.R (person-type "interval" branches).

library(bootnet)

# Avoid creating Rplots.pdf when ggplot_build() / printing touches a device.
pdf(NULL)

set.seed(1)

# 5-variable Gaussian data
p <- 5
n <- 200
data <- as.data.frame(matrix(rnorm(n * p), n, p))
colnames(data) <- paste0("V", seq_len(p))

net <- estimateNetwork(data, default = "pcor")

boot <- suppressMessages(
  bootnet(net, nBoots = 10, type = "person", caseN = 3, verbose = FALSE)
)

## --- plot = "interval" (per-type) ---
g1 <- plot(boot, plot = "interval")
stopifnot(inherits(g1, "ggplot"))
invisible(ggplot2::ggplot_build(g1))

## --- plot = "interval", perNode = TRUE ---
g2 <- plot(boot, plot = "interval", perNode = TRUE)
stopifnot(inherits(g2, "ggplot"))
invisible(ggplot2::ggplot_build(g2))

dev.off()

cat("test-item03-interval-plot.R passed\n")
