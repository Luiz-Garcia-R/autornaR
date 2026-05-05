# R/zzz.R

.onAttach <- function(libname, pkgname) {packageStartupMessage(
    "\n",
    crayon::green("autornaR v.",utils::packageVersion("autornaR")), " loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined RNA-seq data analysis.\n",
    "GitHub: https://github.com/Luiz-Garcia-R/autornaR\n"
  )
}

