summ_logit <- function(model,
                       k = NULL,
                       pct_train = NULL,
                       PI_col = TRUE,
                       or_col = FALSE) {
  
  stopifnot(inherits(model, "glm"), family(model)$family == "binomial")
  set.seed(123)
  
  # -- silent, robust package ensure ------------------------------
  pkgs <- c("jtools","dplyr","tibble","margins","glmtoolbox","finalfit")
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
  
  # --- helpers ---
  inv_logit <- function(eta) 1/(1+exp(-eta))
  tjur_r2   <- function(y, p) { y <- as.numeric(y); mean(p[y==1]) - mean(p[y==0])}
  align_num_string <- function(x, digits = 2) {    # fixed-width numeric, 2dp
    num <- suppressWarnings(as.numeric(x))
    s   <- ifelse(is.na(num), "", sprintf(paste0("%.", digits, "f"), num))
    s   <- ifelse(startsWith(s,"-"), s, paste0("\u2007", s))   # figure-space for positives
    w   <- max(nchar(s, type="width"), na.rm = TRUE)
    paste0(strrep("\u2007", pmax(0L, w - nchar(s, type="width"))), s)
  }
  fw <- function(chars) {                            # fixed-width for p (chars)
    w <- max(nchar(chars, type="width"), na.rm = TRUE)
    paste0(strrep("\u2007", w - nchar(chars, type="width")), chars)
  }
  # --- aGVIF from glmtoolbox on an OLS stub, returned at the COEFFICIENT level ---
  .agvif_glmtoolbox_vec <- function(fitted_model) {
    mf <- stats::model.frame(fitted_model)
    tt <- stats::terms(fitted_model)
    lm_stub <- stats::lm(tt, data = mf)
    
    # suppress gvif()’s printed output
    gvif_out <- suppressMessages(capture.output(G <- glmtoolbox::gvif(lm_stub)))
    
    a_by_term <- as.numeric(G[, 3])^2 # square to put on regular VIF scale
    names(a_by_term) <- rownames(G)
    
    X_full <- stats::model.matrix(lm_stub)
    has_int <- isTRUE(attr(stats::terms(lm_stub), "intercept") == 1)
    coef_all <- names(stats::coef(lm_stub))
    coef_names <- if (has_int) coef_all[-1] else coef_all
    
    assign_full <- attr(X_full, "assign")
    assign_vec  <- if (has_int) assign_full[-1] else assign_full
    term_labels <- attr(stats::terms(lm_stub), "term.labels")
    coef_terms  <- term_labels[assign_vec]
    
    v <- setNames(a_by_term[coef_terms], coef_names)
    if (has_int) v <- c("(Intercept)" = NA_real_, v)
    v
  }
  
  validate_tjur <- function(m, k, pct_train) {
    dat <- stats::model.frame(m)
    y_raw <- model.response(dat)
    n  <- nrow(dat)
    
    # ---- robust 0/1 coding for Y ----
    encode_y01 <- function(y) {
      if (is.factor(y)) {
        lev <- levels(y)
        if (length(lev) != 2L) stop("Outcome must be binary for Tjur's R^2.")
        as.numeric(y == lev[2L])  # second level = 1
      } else if (is.logical(y)) {
        as.numeric(y)
      } else if (is.character(y)) {
        u <- sort(unique(y))
        if (length(u) != 2L) stop("Outcome must be binary for Tjur's R^2.")
        as.numeric(y == u[2L])
      } else {
        as.numeric(y)  # assumes already 0/1 or similar
      }
    }
    
    y <- encode_y01(y_raw)
    
    tjur_fun <- function(y_vec, p_vec) {
      y_vec <- as.numeric(y_vec)
      m1 <- mean(p_vec[y_vec == 1], na.rm = TRUE)
      m0 <- mean(p_vec[y_vec == 0], na.rm = TRUE)
      m1 - m0
    }
    
    # ---------- k-fold CV (aggregated, summ_lm-style) ----------
    if (!is.null(k)) {
      # allow k = "n"
      if (is.character(k)) {
        if (identical(k, "n")) {
          k <- n
        } else {
          stop("`k` must be an integer or 'n'.")
        }
      }
      k <- as.integer(k)
      if (k < 2L || k > n) stop("`k` must be between 2 and n (or 'n').")
      
      folds <- sample(rep_len(seq_len(k), n))
      p_all <- numeric(n)
      
      for (fold in seq_len(k)) {
        tr <- dat[folds != fold, , drop = FALSE]
        te <- dat[folds == fold, , drop = FALSE]
        mi <- stats::glm(stats::formula(m), data = tr, family = binomial())
        p_all[folds == fold] <- stats::predict(mi, newdata = te, type = "response")
      }
      
      val <- tjur_fun(y, p_all)
      return(list(kind = "CV", val = val, k = k, n_test = NA_integer_))
    }
    
    # ---------- single train/test split ----------
    if (!is.null(pct_train) && pct_train > 0 && pct_train < 1) {
      ntr <- floor(pct_train * n)
      idx <- sample.int(n, ntr)
      tr  <- dat[idx, , drop = FALSE]
      te  <- dat[-idx, , drop = FALSE]
      
      mi  <- stats::glm(stats::formula(m), data = tr, family = binomial())
      p_te <- stats::predict(mi, newdata = te, type = "response")
      
      y_te <- encode_y01(model.response(te))
      val  <- tjur_fun(y_te, p_te)
      
      return(list(kind = "Holdout", val = val, k = NA_integer_, n_test = nrow(te)))
    }
    
    NULL
  }
  
  summ_logit_helper <- function(logit_model) {
    # Extract model summary
    data <- logit_model[["model"]]
    summ_logit <- summary(logit_model)
    n <- nrow(data)
    y_var <- logit_model[["formula"]][[2]]
    
    # Number of predictors (excluding intercept)
    p <- length(logit_model[["coefficients"]]) - 1
    
    # Extract parameter names, coefficients, and p-values (excluding intercept)
    parameters_logit <- rownames(summ_logit[["coefficients"]])[-1]
    estimates_logit  <- exp(summ_logit[["coefficients"]][-1, 1])  # Odds Ratio
    p_values_logit   <- summ_logit[["coefficients"]][-1, 4]       # P-values
    intercept_logit  <- exp(summ_logit[["coefficients"]][1, 1])   # Intercept OR
    intercept_p_logit <- summ_logit[["coefficients"]][1, 4]       # Intercept p-value
    
    # Cohen's d (logit-scale -> d)
    cohens_d_logit <- (sqrt(3)/pi) * summ_logit[["coefficients"]][-1, 1]
    cohens_d_logit_intercept <- (sqrt(3)/pi) * summ_logit[["coefficients"]][1, 1]
    
    # Average Marginal Effects
    factor_column <- data.frame(factor = parameters_logit)
    ame_df <- summary(margins(logit_model))[,1:2]
    ame <- dplyr::left_join(factor_column, ame_df, by="factor")$AME
    intercept_probability <- intercept_logit / (1 + intercept_logit)
    
    ## VIF (aGVIF) 
    vifs <- if (ncol(logit_model[["model"]]) > 2) {
      .agvif_glmtoolbox_vec(logit_model)
    } else {
      setNames(rep(NA_real_, length(coef(logit_model))), names(coef(logit_model)))
    }
    
    # Compute Tjur’s R² for Logistic Regression
    p_hat_tjur <- predict(logit_model, type = "response")
    p_hat_1_tjur <- mean(p_hat_tjur[data[[y_var]] == 1])
    p_hat_0_tjur <- mean(p_hat_tjur[data[[y_var]] == 0])
    tjur_r2_logit <- p_hat_1_tjur - p_hat_0_tjur
    
    # Adjusted Tjur’s R² (OLS-style adjustment)
    tjur_r2_adj_logit <- 1 - ((1 - tjur_r2_logit) * (n - 1) / (n - p - 1))
    
    # Format predictor rows
    predictor_rows_logit <- data.frame(
      parameters_logit = parameters_logit,
      estimates_logit  = finalfit::round_tidy(estimates_logit, 2),
      ame_logit        = finalfit::round_tidy(ame, 2),
      p_values_logit   = finalfit::round_tidy(p_values_logit, 3),
      cohens_d_logit   = finalfit::round_tidy(cohens_d_logit, 2),
      vif              = finalfit::round_tidy(vifs[-1], 2)
    )
    
    # Format intercept row
    intercept_row_logit <- data.frame(
      parameters_logit = "(Intercept)",
      estimates_logit  = finalfit::round_tidy(intercept_logit, 2),
      ame_logit        = finalfit::round_tidy(intercept_probability, 2),
      p_values_logit   = finalfit::round_tidy(intercept_p_logit, 3),
      cohens_d_logit   = finalfit::round_tidy(cohens_d_logit_intercept, 2),
      vif              = ""
    )
    
    # R² rows
    r2_row_logit <- data.frame(
      parameters_logit = "Tjur's R²",
      estimates_logit  = finalfit::round_tidy(tjur_r2_logit, 2),
      ame_logit = "", p_values_logit = "", cohens_d_logit = "", vif = ""
    )
    adj_r2_row_logit <- data.frame(
      parameters_logit = "Adjusted Tjur's R²",
      estimates_logit  = finalfit::round_tidy(tjur_r2_adj_logit, 2),
      ame_logit = "", p_values_logit = "", cohens_d_logit = "", vif = ""
    )
    
    # n row
    n_obs_row_logit <- data.frame(
      parameters_logit = "Observations",
      estimates_logit  = n,
      ame_logit = "", p_values_logit = "", cohens_d_logit = "", vif = ""
    )
    
    logit_table <- rbind.data.frame(
      intercept_row_logit,
      predictor_rows_logit,
      r2_row_logit,
      adj_r2_row_logit,
      n_obs_row_logit
    )
    
    colnames_df <- c(
      paste0("DV: ", y_var), 
      "Odds_&_OR", "AME", "p_value", "Cohens_d", "VIF"
    )
    
    rownames(logit_table) <- NULL
    colnames(logit_table) <- colnames_df
    
    logit_table
  }
  
  # --- base objects ---
  mf   <- stats::model.frame(model)
  y    <- as.numeric(model.response(mf))
  X_mm <- stats::model.matrix(model, mf)
  
  # assign/terms for interaction-safe mapping
  asg  <- attr(X_mm, "assign")
  labs <- attr(terms(model), "term.labels")
  col_terms <- ifelse(asg == 0, "(Intercept)", labs[asg])
  names(col_terms) <- colnames(X_mm)
  row_terms <- col_terms[match(names(coef(model)), names(col_terms))]
  names(row_terms) <- names(coef(model))
  
  # --- call helper for OR/AME/d/VIF columns ---
  base_tab   <- summ_logit_helper(model)
  coef_names <- names(coef(model))
  coef_rows  <- match(coef_names, base_tab[[1]])
  stopifnot(all(!is.na(coef_rows)))
  
  OR_col   <- base_tab$`Odds_&_OR`[coef_rows]
  AME_col  <- base_tab$AME[coef_rows]
  p_col    <- base_tab$p_value[coef_rows]
  d_col    <- base_tab$Cohens_d[coef_rows]
  VIF_col  <- base_tab$VIF[coef_rows]
  
  # --- raw Cohen's d for Pratt (numeric, including intercept) ---
  sm_full <- summary(model)$coefficients
  d_raw   <- (sqrt(3)/pi) * sm_full[, "Estimate"]
  names(d_raw) <- rownames(sm_full)
  
  # --- model fit (Tjur’s R²) + optional validation ---
  p_train    <- stats::predict(model, type = "response")
  tjur_train <- tjur_r2(y, p_train)
  vinfo      <- validate_tjur(model, k, pct_train)
  
  # choose target Tjur R² for importance scaling
  if (is.null(vinfo)) {
    R2_target <- tjur_train
    lift_col_name <- "PI"
  } else if (identical(vinfo$kind, "Holdout")) {
    R2_target <- vinfo$val
    lift_col_name <- "PI_Test"
  } else {  # CV
    R2_target <- vinfo$val
    lift_col_name <- "PI_CV"
  }
  
  # --- Cohen's d × raw r(X,Y) Pratt-style importance (term-level) ---
  coef_all <- colnames(X_mm)
  coef_noi <- coef_all[coef_all != "(Intercept)"]
  X_noi    <- X_mm[, coef_noi, drop = FALSE]
  
  r_vec <- vapply(seq_along(coef_noi), function(j) {
    xj <- X_noi[, j]
    if (all(is.na(xj)) || stats::sd(xj) == 0) return(0)
    suppressWarnings(stats::cor(xj, y))
  }, numeric(1))
  names(r_vec) <- coef_noi
  
  d_vec <- d_raw[coef_noi]
  d_vec[is.na(d_vec)] <- 0
  
  raw_contrib <- d_vec * r_vec
  
  mm_assign   <- attr(X_mm, "assign")
  term_labels <- attr(stats::terms(model), "term.labels")
  term_idx    <- mm_assign[coef_all != "(Intercept)"]
  
  pratt_coef_df <- tibble::tibble(
    variable = coef_noi,
    raw      = raw_contrib,
    term     = term_labels[term_idx]
  )
  
  term_sum <- pratt_coef_df |>
    dplyr::group_by(term) |>
    dplyr::summarise(S = sum(raw, na.rm = TRUE), .groups = "drop")
  
  S_tot <- sum(term_sum$S, na.rm = TRUE)
  if (is.finite(S_tot) && S_tot != 0 && is.finite(R2_target) && !is.na(R2_target)) {
    term_sum <- term_sum |>
      dplyr::mutate(Pratt_term = S / S_tot * R2_target)
  } else {
    term_sum$Pratt_term <- NA_real_
  }
  
  term_to_pratt <- setNames(term_sum$Pratt_term, term_sum$term)
  pratt_coef_df$Pratt_coef <- term_to_pratt[pratt_coef_df$term]
  
  pratt_vec <- setNames(rep(NA_real_, length(coef_all)), coef_all)
  pratt_vec[coef_noi] <- pratt_coef_df$Pratt_coef[
    match(coef_noi, pratt_coef_df$variable)
  ]
  
  # --- build a jtools summ object just for printing skeleton ---
  s_obj <- jtools::summ(model, exp = TRUE, confint = FALSE,
                        ci.level = 0.95, digits = 2)
  
  ct <- data.frame(
    AME      = AME_col,
    Cohens_d = d_col,
    VIF      = VIF_col,
    stringsAsFactors = FALSE, row.names = coef_names
  )
  if (isTRUE(or_col)) ct <- cbind(OR = OR_col, ct)
  
  if (isTRUE(PI_col)) {
    PI_vals <- as.numeric(pratt_vec[coef_names])
    PI_vals[1] <- NA_real_  # intercept
    ct[[lift_col_name]] <- PI_vals
  }
  
  # p: force 3 decimals + stars and fixed width
  sm <- coef(summary(model)); p_raw <- as.numeric(sm[coef_names, "Pr(>|z|)"])
  add_stars <- function(p) {
    if (is.na(p)) return("\u2007")
    if (p < 0.001) return("\u2007````")
    if (p < 0.01)  return("\u2007```")
    if (p < 0.05)  return("\u2007``")
    "\u2007"
  }
  
  # explicit "NA" for VIFs that are NA
  if ("VIF" %in% names(ct)) {
    ct$VIF[is.na(ct$VIF)] <- "NA"
  }
  
  # numeric columns fixed-width (including Pratt column)
  pi_cols <- grep("^PI", names(ct), value = TRUE)
  num_set <- intersect(c("OR","AME","Cohens_d", "VIF", pi_cols), names(ct))
  if ("VIF" %in% names(ct) && !is.numeric(ct$VIF)) {
    num_set <- setdiff(num_set, "VIF")
  }
  for (cc in num_set) ct[[cc]] <- align_num_string(ct[[cc]], 2)
  
  # ---- center headers to true visual widths ----
  col_w <- vapply(names(ct), function(nm) {
    vals <- as.character(ct[[nm]])
    vals <- gsub("\u2007", " ", vals, fixed = TRUE)
    max(nchar(c(vals, nm), type = "width"), na.rm = TRUE)
  }, integer(1))
  
  center_header <- function(nm, width_target) {
    cur_w <- nchar(nm, type = "width")
    pad_total <- max(0L, width_target - cur_w)
    left_pad  <- floor(pad_total / 2) + 1
    paste0(strrep("\u2007", left_pad), nm)
  }
  
  names(ct) <- mapply(center_header, names(ct), col_w, USE.NAMES = FALSE)
  
  ct$p_value <- paste0(finalfit::round_tidy(p_raw, 3), sapply(p_raw, add_stars))
  
  s_obj[["coeftable"]] <- as.matrix(ct)
  
  # --- render and edit header blocks like summ_lm() ---
  L <- capture.output(print(s_obj))
  L <- sub("^Dependent Variable:", "Outcome Variable:", L)
  L <- sub("^Type: Generalized linear model", "Type: Logistic Regression", L)
  L <- L[!grepl("^\\s*Family:|^\\s*Link function:", L)]
  L <- L[!grepl("^Standard errors:", L)]
  
  # Replace MODEL FIT section with Tjur's R² (+ validation if given)
  fit_i <- which(grepl("^MODEL FIT:", L))
  if (length(fit_i)==1) {
    after <- fit_i + 1L
    while (after <= length(L) && nzchar(trimws(L[after]))) after <- after + 1L
    r2 <- "\u00B2"
    if (is.null(vinfo)) {
      fit_lines <- c(
        "MODEL FIT:",
        sprintf("Tjur's R%s = %.2f", r2, tjur_train),
        sprintf("Adjusted Tjur's R%s = %.2f", r2,
                1 - ((1 - tjur_train) * (nrow(mf) - 1) /
                       (nrow(mf) - length(coef(model)) - 1)))
      )
    } else if (identical(vinfo$kind, "Holdout")) {
      fit_lines <- c(
        "MODEL FIT:",
        sprintf("Tjur's R%s = %.2f", r2, tjur_train),
        sprintf("Tjur's R%s (test) = %.2f", r2, vinfo$val),
        sprintf("Test set size = %.1f%%", (vinfo$n_test/nrow(mf))*100)
      )
    } else { # CV
      fit_lines <- c(
        "MODEL FIT:",
        sprintf("Tjur's R%s = %.2f", r2, tjur_train),
        sprintf("Tjur's R%s (CV) = %.2f", r2, vinfo$val),
        sprintf("Folds = %s", as.character(vinfo$k))
      )
    }
    L <- c(L[seq_len(fit_i - 1L)], fit_lines, L[seq(from = after, to = length(L))])
  }
  
  ## ---- Force-correct p-value stars post-render (p always "0.xxx") ----
  p_stars <- function(pv_num) {
    if (is.na(pv_num)) return(" ")
    if (pv_num < 0.001) " ***"
    else if (pv_num < 0.01) " **"
    else if (pv_num < 0.05) " *"
    else " "
  }
  
  fix_p_line <- function(line) {
    m <- gregexpr("\\b0\\.\\d{3}\\b(\\s\\*{0,3})?", line, perl = TRUE)
    pos <- m[[1]]
    if (length(pos) == 1 && pos[1] == -1) return(line)
    
    i <- tail(which(pos > 0), 1)
    start <- pos[i]
    len   <- attr(m[[1]], "match.length")[i]
    
    tok <- substr(line, start, start + len - 1)
    num_str <- sub("\\s\\*{0,3}$", "", tok)
    pv <- suppressWarnings(as.numeric(num_str))
    
    stars <- p_stars(pv)
    replacement <- paste0(num_str, stars)
    
    paste0(substr(line, 1, start - 1), replacement,
           substr(line, start + len, nchar(line)))
  }
  
  L <- vapply(L, fix_p_line, character(1))
  L <- gsub("`+", "", L, perl = TRUE)
  
  cat(paste0(paste(L, collapse = "\n"), "\n"))
  cat("p: *** is <.001, ** is <.01, * is <.05\n",
      "|d|: Trivial <.1 < Tiny <.2 < Small <.5 < Medium <.8 < Large <1.4 < Huge\n\n", sep = "")
  
  invisible(s_obj)
}