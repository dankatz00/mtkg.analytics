summ_glht <- function(glht_obj, p_adjust_method = "single-step") {
  
  # -- silent, robust package ensure ------------------------------
  pkgs <- c("dplyr", "Hmisc", "multcomp")
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
  
  round0 <- function(.data, digits = 2) {
    .data[] <- lapply(.data, function(x) {
      if (is.numeric(x)) formatC(round(x, digits), format = "f", digits = digits) else x
    })
    .data
  }
  
  # Get summary with specified adjustment method
  p_method <- tolower(p_adjust_method)
  summ <- summary(glht_obj, test = adjusted(p_method))
  estimates <- summ$test$coefficients
  pvals <- round(summ$test$pvalues, 3)
  
  # Residual SD
  resid_sd <- sigma(glht_obj$model)
  
  # Compute Cohen's d
  cohens_d <- estimates / resid_sd
  
  # Function to add significance stars
  add_stars = function(p) {
    if (is.na(p)) return("")  # Handle missing p-values
    if (p < 0.001) return(" ***")
    if (p < 0.01) return(" **")
    if (p < 0.05) return(" *")
    return("")
  }
  
  # Apply stars to p-values
  p_with_stars = paste0(round0(pvals, 3), sapply(pvals, add_stars))
  
  # Interpret effect size
  label <- case_when(
    abs(cohens_d) < 0.1 ~ "Trivial",
    abs(cohens_d) < 0.2 ~ "Tiny",
    abs(cohens_d) < 0.5 ~ "Small",
    abs(cohens_d) < 0.8 ~ "Medium",
    abs(cohens_d) < 1.4 ~ "Large",
    TRUE ~ "Huge"
  )
  
  # Raw contrast labels
  raw_labels <- names(estimates)
  
  # Get model variables and data
  model_vars <- colnames(glht_obj$model[["model"]])
  model_data <- model.frame(glht_obj$model)
  
  # Identify which variables are factors
  categorical_vars <- model_vars[vapply(model_data[model_vars], function(x) is.factor(x) || is.character(x), logical(1))]
  
  
  # FOR DEBUG
  #print(categorical_vars)
  
  # Build replacements
  clean_labels <- raw_labels
  for (v in categorical_vars) {
    lvls <- levels(as.factor(model_data[[v]]))
    for (lvl in lvls) {
      pattern <- paste0(v, lvl)
      # FOR DEBUG
      #cat("Replacing:", pattern, "→", lvl, "\n")  # Log what should be replaced
      clean_labels <- gsub(pattern, lvl, clean_labels, fixed = TRUE)
    }
  }
  
  # FOR DEBUG
  #print(clean_labels) 
  
  # Final formatting
  clean_labels <- clean_labels %>%
    gsub("-", " - ", .) %>%
    gsub("\\+", " + ", .) %>%
    gsub("=", " = ", .) %>%
    gsub("\\s+", " ", .) %>%
    trimws() %>%
    paste("= 0")
  
  # Determine p-value column label
  p_label <- if (p_method == "none") "p (raw)" else "p (adj)"
  
  # Build data frame
  df <- data.frame(
    Contrast = clean_labels,
    Estimate = round(estimates, 2),
    Cohens_d = round(cohens_d, 2),
    Interpretation = label,
    stringsAsFactors = FALSE
  )
  
  df[[p_label]] <- p_with_stars  # Add p-values as dynamic column name
  
  # Move p column just after Estimate
  df <- df[c("Contrast", "Estimate", p_label, "Cohens_d", "Interpretation")]
  
  # Strip backticks from term names in output
  df$Contrast <- gsub("`", "", df$Contrast)
  
  rownames(df) <- NULL
  
  Hmisc::print.char.matrix(as.matrix(df),
                           col.names = TRUE,
                           row.names = TRUE,
                           hsep = "  ", 
                           vsep = "", 
                           csep = "",
                           col.name.align = "left", 
                           col.txt.align = "left", 
                           cell.align = "left")
  cat("\n")
}