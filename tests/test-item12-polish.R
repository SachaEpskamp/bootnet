# Regression tests for audit item 12: documentation / API-surface polish.
#
# Covered here:
#   - binarize() emits a message (not a warning) when splitting continuous
#     data into binary form.
#   - differenceTest() on an undirected (pcor) bootstrap still resolves edge
#     endpoints via the "--" separator (guards the undirected path after the
#     "--"/"->" both-separator change).
#   - estimateNetwork(default = "Borsboom") now errors via match.arg()
#     instead of returning the removed easter-egg value 42.
#   - Building a case-drop area plot no longer opens a graphics device as a
#     side effect (previously par("pin") created an Rplots.pdf in scripts).

library(bootnet)

set.seed(1)

## ---------------------------------------------------------------------------
## (a) binarize(): median split emits a message, not a warning.
## ---------------------------------------------------------------------------
contData <- as.data.frame(matrix(rnorm(40 * 4), ncol = 4))
colnames(contData) <- paste0("V", 1:4)

msgs <- character(0)
warns <- character(0)
res <- withCallingHandlers(
    bootnet:::binarize(contData, verbose = TRUE),
    message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
    },
    warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)
stopifnot(any(grepl("Splitting data by", msgs, fixed = TRUE)))
# The splitting notice must no longer be a warning:
stopifnot(!any(grepl("Splitting data by", warns, fixed = TRUE)))
# Result is 0/1 binary:
stopifnot(all(unlist(res) %in% c(0, 1)))

## Constant-column guard: a column that is constant after a median split
## should now trigger a real warning.
constData <- data.frame(
    A = c(1, 2, 3, 4, 5, 6),   # normal, splits into 0/1
    B = c(5, 5, 5, 5, 5, 5)    # constant -> constant after any split
)
warns2 <- character(0)
invisible(withCallingHandlers(
    bootnet:::binarize(constData, verbose = FALSE),
    warning = function(w) {
        warns2 <<- c(warns2, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
))
stopifnot(any(grepl("constant", warns2, ignore.case = TRUE)))

## ---------------------------------------------------------------------------
## (b) differenceTest() on an undirected pcor bootstrap ("--" IDs).
## ---------------------------------------------------------------------------
set.seed(2)
gdata <- as.data.frame(matrix(rnorm(60 * 4), ncol = 4))
colnames(gdata) <- paste0("V", 1:4)

boot <- bootnet(gdata, default = "pcor", nBoots = 8, verbose = FALSE,
                statistics = c("edge", "strength"))

# Strength difference test (node IDs, no separator involved):
dtStrength <- differenceTest(boot, x = "V1", y = "V2",
                             measure = "strength", verbose = FALSE)
stopifnot(is.data.frame(dtStrength))
stopifnot(nrow(dtStrength) >= 1)

# Edge difference test using x2/y2 endpoints -> exercises the "--" branch of
# the both-separator resolution added in item 12.
dtEdge <- differenceTest(boot, x = "V1", x2 = "V2", y = "V3", y2 = "V4",
                         measure = "edge", verbose = FALSE)
stopifnot(is.data.frame(dtEdge))
stopifnot(nrow(dtEdge) >= 1)
# Directed "->" path is not exercised here: the cheap defaults (pcor) produce
# undirected networks only; directed defaults (e.g. relimp) are heavy, so the
# directed separator is left to the code path guarded by "--" symmetry.

## ---------------------------------------------------------------------------
## (c) estimateNetwork(default = "Borsboom") now errors (no easter egg).
## ---------------------------------------------------------------------------
borsboom <- tryCatch(
    estimateNetwork(gdata, default = "Borsboom", verbose = FALSE),
    error = function(e) e
)
stopifnot(inherits(borsboom, "error"))
stopifnot(!identical(borsboom, 42))

## ---------------------------------------------------------------------------
## (d) No-device safety: building a case-drop area plot must not create an
##     Rplots.pdf via par("pin").
## ---------------------------------------------------------------------------
owd <- setwd(tempdir())
on.exit(setwd(owd), add = TRUE)

# Make sure we start clean and that no device is currently open.
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
while (!is.null(grDevices::dev.list())) grDevices::dev.off()

set.seed(3)
bootPerson <- bootnet(gdata, default = "pcor", type = "case",
                      nBoots = 8, caseN = 3, verbose = FALSE,
                      statistics = "strength")

# Default plot for a case-drop bootstrap is the "area" branch that used
# par("pin"). ggplot_build forces construction without opening a device.
g <- plot(bootPerson, statistics = "strength")
invisible(ggplot2::ggplot_build(g))

stopifnot(!file.exists("Rplots.pdf"))

setwd(owd)

cat("item12-polish: OK\n")
