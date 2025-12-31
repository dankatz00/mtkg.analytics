.onLoad <- function(libname, pkgname) {
  
  # allow you to disable auto-install if you ever need to
  if (isFALSE(getOption("mtkg_analytics.auto_install", TRUE))) return(invisible(NULL))
  
  pkgs <- c(
    "dplyr", "tibble", "magrittr",
    "jtools", "finalfit", "glmtoolbox",
    "estimatr", "margins",
    "multcomp", "Hmisc",
    "psych", "lsr", "pwr",
    "survival"
  )
  
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    
    repos <- getOption("repos")
    if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
      repos <- c(CRAN = "https://cloud.r-project.org")
    }
    
    # quiet install
    try(
      suppressWarnings(suppressMessages(
        invisible(capture.output(
          utils::install.packages(missing, quiet = TRUE, repos = repos),
          type = "output"
        ))
      )),
      silent = TRUE
    )
  }
  
  # after install attempt, verify
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    stop(
      "mtkg_analytics could not install/load required packages: ",
      paste(still_missing, collapse = ", ")
    )
  }
  
  # quiet attach so unqualified calls like case_when() and %>% work
  for (p in pkgs) {
    try(suppressWarnings(suppressMessages(
      library(p, character.only = TRUE)
    )), silent = TRUE)
  }
  
  invisible(NULL)
}