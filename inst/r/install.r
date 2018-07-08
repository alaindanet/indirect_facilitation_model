#!/usr/bin/env r

if (is.null(argv) | length(argv)<1) {
  cat("Usage: installr.r pkg1 [pkg2 pkg3 ...]\n")
  q()
}

## adjust as necessary, see help('download.packages')
repos <- "https://cran.rstudio.com"

## this makes sense on Debian where no packages touch /usr/local
lib.loc <- "/home/alain/R/x86_64-pc-linux-gnu-library/3.4"

install.packages(argv, lib.loc, repos)
