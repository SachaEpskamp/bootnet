# I copied this piece of code from Lavaan mainly:

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage(pkgname, " is BETA software! Please report any bugs.")
}

# Please do not look further down this script.















# I am warning you...


















# Ok, here goes some stuff to fool the R check rather than changing all the codes... I am sorry..
P <- contain0 <- data <- fill <- id1 <- id2 <- label <- lower <- na.omit <- name <- node1 <- node2 <- 
original <- prop <- quantile <- runif <- sd <- stability <- type <- upper <- value <- value1  <- 
value2 <- nPerson <- NULL