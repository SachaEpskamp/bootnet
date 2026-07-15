# I copied this piece of code from Lavaan mainly:

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  # packageStartupMessage(pkgname, " is BETA software! Please report any bugs.")
  packageStartupMessage("For questions and issues, please see github.com/SachaEpskamp/bootnet.")
}
