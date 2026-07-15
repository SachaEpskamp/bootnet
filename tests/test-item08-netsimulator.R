# Regression test for audit item 8: cleanup of the vendored parSim engine.
#
# Fixes in R/parSim.R:
#   - dropped top-level library() calls (rely on NAMESPACE imports)
#   - renamed parSim -> bootnet_parSim (no longer shadows the parSim package)
#   - fixed the inverted 'exclude' filter: matching conditions are now EXCLUDED
#
# This exercises netSimulator() end to end and tests the exclude fix directly on
# the (internal) bootnet_parSim, since netSimulator()'s public API does not
# expose 'exclude'.

library(bootnet)

# The internal engine must not leak into the global namespace and must not
# shadow the parSim package.
stopifnot(!exists("parSim"))
stopifnot(exists("bootnet_parSim", where = asNamespace("bootnet"), inherits = FALSE))

set.seed(1)

# --- netSimulator smoke test -------------------------------------------------
sim <- netSimulator(
  input = genGGM(4),
  nCases = 50,
  nReps = 2,
  default = "pcor",
  verbose = FALSE
)

stopifnot(is.data.frame(sim))
stopifnot(all(c("sensitivity", "specificity") %in% names(sim)))
# One condition (single nCases value), nReps = 2 -> 2 rows.
stopifnot(nrow(sim) == 2)

# --- exclude fix: matching conditions must be dropped ------------------------
# Trivial expression that just returns the condition value as a data frame.
parSimFun <- bootnet:::bootnet_parSim

set.seed(2)
res <- parSimFun(
  nCases = c(50, 100),
  reps = 2,
  write = FALSE,
  nCores = 1,
  exclude = list(~ nCases == 50),
  expression = data.frame(value = nCases)
)

stopifnot(is.data.frame(res))
# The excluded condition (nCases == 50) must be absent; only nCases == 100 rows
# remain (2 reps).
stopifnot(all(res$nCases == 100))
stopifnot(!any(res$nCases == 50))
stopifnot(nrow(res) == 2)

# Sanity check without exclude: both conditions present.
set.seed(3)
resAll <- parSimFun(
  nCases = c(50, 100),
  reps = 2,
  write = FALSE,
  nCores = 1,
  expression = data.frame(value = nCases)
)
stopifnot(setequal(unique(resAll$nCases), c(50, 100)))
stopifnot(nrow(resAll) == 4)

cat("item08-netsimulator: OK\n")
