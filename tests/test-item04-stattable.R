# Regression tests for audit item 4 (R/bootTable.R, statTable):
#
# (a) "bridgeInDegree" and "bridgeOutDegree" were listed as valid statistics
#     and networktools::bridge() output was renamed to those keys, but there
#     were no storage blocks for them, so requesting them silently produced
#     no rows in the boot/sample tables.
# (b) Copy-paste bug in the hybrid try-error fallback: the NA fallback tibble
#     was assigned to tables$rspbc instead of tables$hybrid, so when both
#     "rspbc" and "hybrid" were requested and hybrid failed, the valid rspbc
#     rows were silently overwritten by NA hybrid rows.

library(bootnet)

## ---------------------------------------------------------------------------
## (a) bridgeInDegree / bridgeOutDegree are stored
## ---------------------------------------------------------------------------
# Note: networktools::bridge() only returns "Bridge Indegree"/"Bridge Outdegree"
# for directed networks, so we force directed = TRUE via bridgeArgs (pcor
# networks are symmetric, in which case bridge() auto-detects undirected and
# omits these measures).

set.seed(42)
n <- 200
data <- as.data.frame(matrix(rnorm(n * 6), nrow = n, ncol = 6))
colnames(data) <- paste0("V", 1:6)
comm <- c("A", "A", "A", "B", "B", "B")

net <- estimateNetwork(data, default = "pcor", verbose = FALSE)

boot <- bootnet(net, nBoots = 2,
                communities = comm,
                statistics = c("bridgeInDegree", "bridgeOutDegree"),
                bridgeArgs = list(directed = TRUE),
                verbose = FALSE)

sampTab <- boot$sampleTable

stopifnot(is.data.frame(sampTab))
stopifnot("bridgeInDegree" %in% sampTab$type)
stopifnot("bridgeOutDegree" %in% sampTab$type)

inRows  <- sampTab[sampTab$type == "bridgeInDegree", ]
outRows <- sampTab[sampTab$type == "bridgeOutDegree", ]

# One row per node, with node labels and non-NA numeric values:
stopifnot(nrow(inRows) == 6, nrow(outRows) == 6)
stopifnot(setequal(inRows$node1, colnames(data)))
stopifnot(setequal(outRows$node1, colnames(data)))
stopifnot(!any(is.na(inRows$value)), !any(is.na(outRows$value)))
stopifnot(is.numeric(inRows$value), is.numeric(outRows$value))

# Bootstrap table must contain them too:
stopifnot("bridgeInDegree" %in% boot$bootTable$type)
stopifnot("bridgeOutDegree" %in% boot$bootTable$type)
stopifnot(!any(is.na(boot$bootTable$value[
  boot$bootTable$type %in% c("bridgeInDegree", "bridgeOutDegree")])))

cat("bridgeInDegree/bridgeOutDegree storage test passed\n")

## ---------------------------------------------------------------------------
## (b) hybrid try-error fallback must not clobber rspbc rows
## ---------------------------------------------------------------------------
# statTable() is internal; call it directly the way bootnet() does
# (see R/bootnet.R: statTable(sampleResult, name = "sample", ...)).
# NetworkToolbox::hybrid() does not reliably error on small graphs (and on
# graphs where it does error, rspbc errors too, since hybrid(BC = "random")
# calls rspbc internally), so we deterministically force the try-error branch
# by temporarily masking hybrid() in the NetworkToolbox namespace.

set.seed(7)
data2 <- as.data.frame(matrix(rnorm(200 * 3), nrow = 200, ncol = 3))
colnames(data2) <- paste0("V", 1:3)
net2 <- estimateNetwork(data2, default = "pcor", verbose = FALSE)
stopifnot(is(net2, "bootnetResult"))

ns <- asNamespace("NetworkToolbox")
orig_hybrid <- get("hybrid", envir = ns)
unlockBinding("hybrid", ns)
assign("hybrid", function(...) stop("forced hybrid failure for test"), envir = ns)

tab <- tryCatch(
  bootnet:::statTable(net2, name = "sample",
                      statistics = c("rspbc", "hybrid")),
  finally = {
    assign("hybrid", orig_hybrid, envir = ns)
    lockBinding("hybrid", ns)
  }
)

rspbcRows  <- tab[tab$type == "rspbc", ]
hybridRows <- tab[tab$type == "hybrid", ]

# rspbc rows must survive with non-NA values:
stopifnot(nrow(rspbcRows) == 3)
stopifnot(!any(is.na(rspbcRows$value)))

# hybrid rows must be present as the NA fallback:
stopifnot(nrow(hybridRows) == 3)
stopifnot(all(is.na(hybridRows$value)))

# Sanity: hybrid mask has been restored:
stopifnot(identical(get("hybrid", envir = ns), orig_hybrid))

cat("hybrid/rspbc fallback test passed\n")

cat("test-item04-stattable.R passed\n")
