cohens_effsize <- function(test, print_test = TRUE) {
  # -- silent, robust package ensure --------------------------------
  pkgs <- c("lsr", "pwr")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      repos <- getOption("repos")
      if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
        repos <- c(CRAN = "https://cloud.r-project.org")
      }
      try(suppressWarnings(suppressMessages(
        invisible(capture.output(utils::install.packages(pkg, quiet = TRUE, repos = repos), type = "output"))
      )), silent = TRUE)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' could not be loaded even after silent install.")
    }
  }
  
  # -- helpers ----------------------------------------------------------------
  .interp_word <- function(v) {
    a <- abs(v)
    if (a < 0.1) "Trivial"
    else if (a < 0.2) "Tiny"
    else if (a < 0.5) "Small"
    else if (a < 0.8) "Medium"
    else if (a < 1.4) "Large"
    else "Huge"
  }
  .get_obj <- function(name_str) {
    nm  <- trimws(name_str)
    obj <- tryCatch(eval(str2lang(nm), envir = parent.frame()),
                    error = function(e) get(nm, envir = parent.frame()))
    if (is.data.frame(obj) || is.matrix(obj)) obj <- obj[, 1, drop = TRUE]
    if (is.list(obj) && length(obj) == 1) obj <- obj[[1]]
    obj <- unlist(obj, use.names = FALSE)
    if (is.factor(obj)) obj <- as.numeric(as.character(obj))
    if (is.character(obj)) obj <- suppressWarnings(as.numeric(obj))
    obj
  }
  .get_prop <- function(est, key_primary, key_alt = NULL) {
    if (!is.null(est[[key_primary]])) return(unname(est[[key_primary]]))
    if (!is.null(key_alt) && !is.null(est[[key_alt]])) return(unname(est[[key_alt]]))
    stop("Expected proportion '", key_primary,
         if (!is.null(key_alt)) paste0("' or '", key_alt, "'") else "", " not found.")
  }
  .print_result <- function(label, value) {
    value <- round(value, 2)
    cat(sprintf("%s = %.2f (%s Effect)\n", label, value, .interp_word(value)))
  }
  .ordinal <- function(n) {
    # n is integer 0-100
    if (n %% 100L %in% c(11L, 12L, 13L)) suf <- "th"
    else if (n %% 10L == 1L) suf <- "st"
    else if (n %% 10L == 2L) suf <- "nd"
    else if (n %% 10L == 3L) suf <- "rd"
    else suf <- "th"
    paste0(n, suf)
  }
  
  # -- main -------------------------------------------------------------------
  method   <- test[["method"]]
  data_str <- test[["data.name"]]
  
  if (print_test) { print(test) }
  
  if (identical(method, "Welch Two Sample t-test")) {
    # A vs B (order taken from data.name)
    parts <- strsplit(data_str, " and ", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Could not parse data.name into two objects: ", data_str)
    x <- .get_obj(parts[1]); y <- .get_obj(parts[2])
    
    d_signed <- lsr::cohensD(x = x, y = y, method = "unequal")
    .print_result("Cohen's d", d_signed)
    
  } else if (identical(method, "One Sample t-test")) {
    x  <- .get_obj(data_str)
    mu <- unname(test[["null.value"]][["mean"]])
    
    d_signed <- lsr::cohensD(x = x, mu = mu)
    .print_result("Cohen's d", d_signed)
    
  } else if (identical(method, "2-sample test for equality of proportions with continuity correction")) {
    est <- test[["estimate"]]
    p1 <- .get_prop(est, "prop 1", "p1")
    p2 <- .get_prop(est, "prop 2", "p2")
    
    h_signed <- pwr::ES.h(p1 = p1, p2 = p2)    # keep sign for percentile
    .print_result("Cohen's h", abs(h_signed))
    
  } else if (identical(method, "1-sample proportions test with continuity correction")) {
    est_p  <- unname(test[["estimate"]][["p"]])
    null_p <- unname(test[["null.value"]][["p"]])
    
    h_signed <- pwr::ES.h(p1 = est_p, p2 = null_p)
    .print_result("Cohen's h", abs(h_signed))
    
  } else {
    stop("Unsupported test method: ", method)
  }
}