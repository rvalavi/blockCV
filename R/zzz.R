.onAttach <- function(libname, pkgname) {
  packageStartupMessage("blockCV ", utils::packageVersion("blockCV"))
}
