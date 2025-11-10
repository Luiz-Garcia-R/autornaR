# R/zzz.R

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    crayon::green("autornaR "), "loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined RNA-seq data analysis.\n",
    "Use ", crayon::green("?autornaR"), " for general help.\n",
    "--------------------------------------------------\n"
  )
}
