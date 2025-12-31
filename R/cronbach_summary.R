cronbach_summary <- function(alpha_obj) {
  # Make sure psych is installed
  if (!requireNamespace("psych", quietly = TRUE)) {
    install.packages("psych")
  }
  # Make sure psych is loaded
  if (!"package:psych" %in% search()) {
    library(psych)
  }
  # Basic check that this looks like a psych::alpha object
  if (is.null(alpha_obj$alpha.drop) || is.null(alpha_obj$total)) {
    stop("`alpha_obj` must be an object returned by psych::alpha().")
  }
  
  # Overall alpha
  overall <- round(alpha_obj$total$raw_alpha, 2)
  
  # Alpha if each item is dropped
  item_alpha <- alpha_obj$alpha.drop[, "raw_alpha"]
  names(item_alpha) <- rownames(alpha_obj$alpha.drop)
  item_alpha <- round(item_alpha, 2)
  
  ## ---- printing (no bracketed indices) ----
  cat("\n-----------------\n")
  cat("Cronbach's alpha:\n")
  cat(overall, "\n\n", sep = "")
  
  cat("Alpha if each item is dropped:\n")
  for (nm in names(item_alpha)) {
    cat(nm, ": ", item_alpha[[nm]], "\n", sep = "")
  }
}