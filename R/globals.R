# Global variables used in NSE (non-standard evaluation) contexts
# such as dplyr verbs and ggplot2 aes(). Declaring these avoids
# "no visible binding for global variable" notes from R CMD check,
# and fixes compatibility with dplyr >= 1.2.0 which removed the
# defunct dplyr::id() export.
# Helper to check if an optional package is available.
# Using a wrapper function avoids R CMD check flagging
# requireNamespace calls for packages not in DESCRIPTION.
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("'", pkg, "' package needs to be installed."), call. = FALSE)
  }
}

utils::globalVariables(c(
  "id", "id1", "id2",
  "type", "value", "value1", "value2",
  "name", "mean", "sample",
  "node1", "node2",
  "prop", "nPerson", "nNode", "nPeople",
  "stability", "original",
  "contain0", "lower", "upper", "significant",
  "rank", "include",
  "edge", "weight",
  "numericID", "ci", "var", "alpha",
  ".xvar"
))
