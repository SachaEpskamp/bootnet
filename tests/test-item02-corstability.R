# Audit item 2: corStability() field-name mismatch (nNodes -> nNode)
# See R/Smetric.R (corStability), R/bootnet.R (manual parametric branch),
# R/bootTable.R (statTable reads x$nNode).

library(bootnet)

set.seed(1)

## --- Node-drop bootstrap: corStability must return a named numeric vector ---

# 5-variable Gaussian data
p <- 5
n <- 200
data <- as.data.frame(matrix(rnorm(n * p), n, p))
colnames(data) <- paste0("V", seq_len(p))

net <- estimateNetwork(data, default = "pcor")

boot <- suppressMessages(
  bootnet(net, nBoots = 6, type = "node", verbose = FALSE)
)

cs <- suppressMessages(corStability(boot, verbose = FALSE))

stopifnot(is.numeric(cs))
stopifnot(!is.null(names(cs)))
stopifnot(length(cs) >= 1)
stopifnot(all(!is.na(names(cs))))


## --- Manual parametric bootstrap: bootTable$nNode must be populated ---

g <- matrix(0, 4, 4)
g[1, 2] <- g[2, 1] <- 0.3
g[3, 4] <- g[4, 3] <- 0.3
colnames(g) <- rownames(g) <- paste0("V", seq_len(ncol(g)))

pboot <- suppressMessages(
  bootnet(graph = g, sampleSize = 100, intercepts = rep(0, ncol(g)),
          type = "parametric", model = "GGM", default = "pcor",
          nBoots = 3, verbose = FALSE)
)

stopifnot("nNode" %in% colnames(pboot$bootTable))
stopifnot(!is.null(pboot$bootTable$nNode))
stopifnot(all(!is.na(pboot$bootTable$nNode)))
stopifnot(all(pboot$bootTable$nNode == ncol(g)))

# The sampleResult field is named nNode (not nNodes)
stopifnot(!is.null(pboot$sample$nNode))
stopifnot(is.null(pboot$sample$nNodes))

cat("test-item02-corstability.R passed\n")
