summ_lm_robust <- function(model, k = NULL, pct_train = NULL,
                           cohens_d_col = TRUE,
                           sd_d        = FALSE,
                           vif_col     = TRUE,
                           PI_col      = TRUE) {
  
  if (!inherits(model, "lm_robust")) stop("`model` must be a fitted estimatr::lm_robust object.")
  set.seed(123)
  
  # -- silent, robust package ensure ------------------------------
  pkgs <- c("dplyr","glmtoolbox","finalfit","jtools","tibble", "estimatr")
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
  
  ## ---------------- helpers ----------------
  add_stars <- function(p) {
    if (is.na(p)) return("\u2007")
    if (p < 0.001) return("\u2007````")
    if (p < 0.01)  return("\u2007```")
    if (p < 0.05)  return("\u2007``")
    "\u2007"
  }
  
  align_num_string <- function(x) {
    out <- as.character(x)
    was_neg_zero <- grepl("^\\s*-0+(\\.0+)?$", out)
    num <- suppressWarnings(as.numeric(out))
    is_num <- !is.na(num)
    out[is_num] <- sprintf("%.2f", num[is_num])
    out[was_neg_zero] <- "-0.00"
    pos <- is_num & !grepl("^\\s*-", out)
    out[pos] <- paste0("\u2007", out[pos])
    out
  }
  
  calc_R2 <- function(y, yhat, has_int = TRUE, y_ref_mean = NULL) {
    SST <- if (has_int) {
      sum((y - if (is.null(y_ref_mean)) mean(y) else y_ref_mean)^2)
    } else sum(y^2)
    1 - sum((y - yhat)^2) / SST
  }
  
  # Within R2 computed by demeaning within FE group
  calc_within_R2 <- function(y, yhat, g) {
    g <- as.character(g)
    y_tilde    <- y    - ave(y,    g, FUN = mean)
    yhat_tilde <- yhat - ave(yhat, g, FUN = mean)
    1 - sum((y_tilde - yhat_tilde)^2) / sum((y_tilde - mean(y_tilde))^2)
  }
  
  cohen_d_from_t <- function(tval, df, beta_sign) {
    t2  <- tval^2
    pr2 <- t2 / pmax(t2 + df, .Machine$double.eps)
    r   <- sqrt(pmax(pr2, 0))
    d   <- 2 * r / sqrt(pmax(1 - r^2, .Machine$double.eps))
    d * beta_sign
  }
  
  .get_data_from_call <- function(m) {
    if (is.null(m$call$data)) stop("lm_robust object has no `call$data`; cannot reconstruct data.")
    eval(m$call$data, envir = parent.frame())
  }
  
  .strip_tilde <- function(x) sub("^~\\s*", "", x)
  
  .extract_labels <- function(m) {
    fe_label <- NULL
    if (isTRUE(m$fes) && !is.null(m$felevels)) {
      fe_names <- names(m$felevels)
      if (length(fe_names)) fe_label <- paste(fe_names, collapse = ", ")
    }
    
    cl_label <- NULL
    if (isTRUE(m$clustered) && !is.null(m$call$clusters)) {
      cl_expr <- paste(deparse(m$call$clusters, width.cutoff = 500), collapse = "")
      cl_expr <- gsub("\\s+", " ", cl_expr)
      cl_expr <- .strip_tilde(cl_expr)
      
      if (!is.null(m$nclusters) && is.finite(m$nclusters)) {
        cl_label <- paste0(cl_expr, " (", as.integer(m$nclusters), " groups)")
      } else {
        cl_label <- cl_expr
      }
    }
    
    list(fixed_effects = fe_label, clusters = cl_label)
  }
  
  .within_r2 <- function(m) {
    if (isTRUE(m$fes) && !is.null(m$proj_r.squared) && is.finite(m$proj_r.squared)) {
      return(as.numeric(m$proj_r.squared))
    }
    NA_real_
  }
  
  # OLS stub for summ() formatting and VIF (terms only; FE not in terms)
  .ols_stub_from_terms <- function(fml, data_used) {
    stats::lm(formula = fml, data = data_used)
  }
  
  .agvif_glmtoolbox_vec_lm <- function(fitted_lm) {
    G <- {
      suppressMessages(suppressWarnings(
        capture.output(out <- glmtoolbox::gvif(fitted_lm))
      ))
      out
    }
    
    a_by_term <- as.numeric(G[, 3])^2 # square to put on regular VIF scale
    names(a_by_term) <- rownames(G)
    
    X_full   <- stats::model.matrix(fitted_lm)
    has_int  <- isTRUE(attr(stats::terms(fitted_lm), "intercept") == 1)
    coef_all <- names(stats::coef(fitted_lm))
    coef_names <- if (has_int) coef_all[-1] else coef_all
    
    assign_full <- attr(X_full, "assign")
    assign_vec  <- if (has_int) assign_full[-1] else assign_full
    term_labels <- attr(stats::terms(fitted_lm), "term.labels")
    coef_terms  <- term_labels[assign_vec]
    
    v <- setNames(a_by_term[coef_terms], coef_names)
    if (has_int) v <- c("(Intercept)" = NA_real_, v)
    v
  }
  
  ## ---------------- pull specs ----------------
  tt_eval   <- model$terms
  fml       <- stats::formula(tt_eval)
  
  fe_spec   <- if (!is.null(model$call$fixed_effects)) model$call$fixed_effects else NULL
  clus_spec <- if (!is.null(model$call$clusters))      model$call$clusters      else NULL
  se_type   <- if (!is.null(model$se_type)) model$se_type else "HC2"
  
  ## ---------------- build data_used with FE/cluster vars included ----------------
  data_obj <- .get_data_from_call(model)
  
  vars_needed <- unique(c(
    all.vars(fml),
    if (!is.null(fe_spec))   all.vars(fe_spec)   else character(0),
    if (!is.null(clus_spec)) all.vars(clus_spec) else character(0)
  ))
  if (!length(vars_needed)) stop("Could not determine variables needed for lm_robust refits.")
  
  tmp <- data_obj[, vars_needed, drop = FALSE]
  keep <- stats::complete.cases(tmp)
  data_used <- data_obj[keep, , drop = FALSE]
  n <- nrow(data_used)
  if (n < 2) stop("After dropping missing rows, not enough observations remain.")
  
  ## ---------------- evaluated terms-only frame and matrix ----------------
  mf_terms <- stats::model.frame(tt_eval, data = data_used)
  y <- stats::model.response(mf_terms)
  has_intercept <- isTRUE(attr(tt_eval, "intercept") == 1)
  X_full <- stats::model.matrix(tt_eval, data = mf_terms)
  
  ## ---------------- model fit numbers from lm_robust object ----------------
  r2_in  <- if (!is.null(model$r.squared) && is.finite(model$r.squared)) as.numeric(model$r.squared) else NA_real_
  adj_in <- if (!is.null(model$adj.r.squared) && is.finite(model$adj.r.squared)) as.numeric(model$adj.r.squared) else NA_real_
  within_r2 <- .within_r2(model)
  
  ## ---------------- FE variable name (first FE only, matches your examples) ----------------
  fe_group_var <- NULL
  if (isTRUE(model$fes) && !is.null(model$felevels)) {
    fe_names <- names(model$felevels)
    if (length(fe_names) >= 1) fe_group_var <- fe_names[1]
  }
  
  ## ---------------- refit helpers ----------------
  .fit_fold <- function(train_df) {
    args <- list(formula = fml, data = train_df, se_type = se_type)
    
    if (!is.null(clus_spec)) {
      cl_eval <- clus_spec
      if (is.language(cl_eval) || is.name(cl_eval)) {
        cl_eval <- eval(cl_eval, envir = train_df, enclos = parent.frame())
      }
      args$clusters <- cl_eval
    }
    
    if (!is.null(fe_spec)) args$fixed_effects <- fe_spec
    
    do.call(estimatr::lm_robust, args)
  }
  
  .pred_fold <- function(fit, test_df) {
    out <- try(stats::predict(fit, newdata = test_df), silent = TRUE)
    if (inherits(out, "try-error")) stop("predict() failed inside validation.")
    as.numeric(out)
  }
  
  ## ---------------- CV / test R2 (FE => compute WITHIN validation R2) ----------------
  r2_oos <- NA_real_
  fold_id <- NULL
  
  if (!is.null(k)) {
    if (identical(k, "n")) k <- n
    if (k < 2 || k > n) stop("`k` must be an integer between 2 and n.")
    
    if (!is.null(fe_group_var)) {
      g <- as.character(data_used[[fe_group_var]])
      min_sz <- min(table(g))
      if (min_sz < 2) stop("At least one FE group has only 1 observation; CV with FE is not supported.")
      
      fold_id <- integer(n)
      for (lev in unique(g)) {
        idx <- which(g == lev)
        fold_id[idx] <- sample(rep(1:k, length.out = length(idx)))
      }
    } else {
      fold_id <- sample(rep(1:k, length.out = n))
    }
    
    y_pred_all <- numeric(n)
    for (fold in seq_len(k)) {
      tr <- data_used[fold_id != fold, , drop = FALSE]
      te <- data_used[fold_id == fold, , drop = FALSE]
      fit <- .fit_fold(tr)
      y_pred_all[fold_id == fold] <- .pred_fold(fit, te)
    }
    
    if (!is.null(fe_group_var)) {
      g_all <- as.character(data_used[[fe_group_var]])
      r2_oos <- calc_within_R2(y, y_pred_all, g_all)   # CV Within R2
    } else {
      r2_oos <- calc_R2(y, y_pred_all, has_intercept)
    }
  }
  
  if (!is.null(pct_train)) {
    if (pct_train <= 0 || pct_train >= 1) stop("`pct_train` must be in (0, 1).")
    
    idx_tr <- rep(FALSE, n)
    
    if (!is.null(fe_group_var)) {
      g <- as.character(data_used[[fe_group_var]])
      min_sz <- min(table(g))
      if (min_sz < 2) stop("At least one FE group has only 1 observation; train/test with FE is not supported.")
      
      for (lev in unique(g)) {
        idx <- which(g == lev)
        n_tr <- max(1, floor(pct_train * length(idx)))
        idx_tr[sample(idx, size = n_tr)] <- TRUE
      }
    } else {
      idx_tr[sample(seq_len(n), size = floor(pct_train * n))] <- TRUE
    }
    
    tr <- data_used[idx_tr, , drop = FALSE]
    te <- data_used[!idx_tr, , drop = FALSE]
    
    fit_tr <- .fit_fold(tr)
    y_hat_te <- .pred_fold(fit_tr, te)
    
    y_te <- stats::model.response(stats::model.frame(tt_eval, data = te))
    y_tr <- stats::model.response(stats::model.frame(tt_eval, data = tr))
    
    if (!is.null(fe_group_var)) {
      g_te <- as.character(te[[fe_group_var]])
      r2_oos <- calc_within_R2(y_te, y_hat_te, g_te)   # Test Within R2
    } else {
      r2_oos <- calc_R2(y_te, y_hat_te, has_intercept, mean(y_tr))
    }
  }
  
  ## ---------------- OLS stub for summ() + VIF ----------------
  stub_lm <- .ols_stub_from_terms(fml, data_used)
  
  vif_vec <- setNames(rep(NA_real_, length(model$coefficients)), names(model$coefficients))
  if (isTRUE(vif_col)) {
    Xs <- stats::model.matrix(stub_lm)
    has_int_stub <- isTRUE(attr(stats::terms(stub_lm), "intercept") == 1)
    p_stub <- ncol(Xs) - if (has_int_stub) 1L else 0L
    
    if (p_stub >= 2L) {
      v_all <- .agvif_glmtoolbox_vec_lm(stub_lm)
      keep <- names(model$coefficients)
      vif_vec <- v_all[keep]
    }
  }
  
  ## ---------------- robust coefficients + Cohen's d ----------------
  b_hat  <- model$coefficients
  t_hat  <- model$statistic
  df_hat <- model$df
  p_hat  <- model$p.value
  
  beta_sign <- sign(as.numeric(b_hat))
  d_hat <- cohen_d_from_t(as.numeric(t_hat), as.numeric(df_hat), beta_sign)
  names(d_hat) <- names(b_hat)
  
  ## ---------------- PI (Pratt-style; FE => sum to WITHIN (validation if present)) ----------------
  if (isTRUE(model$fes)) {
    if (!is.null(k) || !is.null(pct_train)) {
      R2_target <- if (is.finite(r2_oos)) r2_oos else within_r2
    } else {
      R2_target <- within_r2
    }
  } else {
    if (!is.null(k) || !is.null(pct_train)) {
      R2_target <- if (is.finite(r2_oos)) r2_oos else r2_in
    } else {
      R2_target <- r2_in
    }
  }
  
  coef_all <- colnames(X_full)
  coef_noi <- coef_all[coef_all != "(Intercept)"]
  X_noi <- X_full[, coef_noi, drop = FALSE]
  
  r_vec <- vapply(seq_along(coef_noi), function(j) {
    xj <- X_noi[, j]
    if (all(is.na(xj)) || stats::sd(xj) == 0) return(0)
    suppressWarnings(stats::cor(xj, y))
  }, numeric(1))
  names(r_vec) <- coef_noi
  
  d_vec <- d_hat[match(coef_noi, names(d_hat))]
  d_vec[is.na(d_vec)] <- 0
  raw_contrib <- d_vec * r_vec
  
  pratt_coef_df <- tibble::tibble(variable = coef_noi, raw = raw_contrib)
  
  mm_assign <- attr(X_full, "assign")
  mm_terms  <- attr(tt_eval, "term.labels")
  term_idx  <- mm_assign[coef_all != "(Intercept)"]
  pratt_coef_df$term <- mm_terms[term_idx]
  
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
  
  pi_by_coef <- setNames(rep(NA_real_, length(names(b_hat))), names(b_hat))
  col_to_term <- setNames(mm_terms[term_idx], coef_noi)
  
  for (nm in names(pi_by_coef)) {
    if (nm %in% names(col_to_term)) {
      trm <- col_to_term[[nm]]
      pi_by_coef[[nm]] <- term_to_pratt[[trm]]
    } else if (nm %in% names(term_to_pratt)) {
      pi_by_coef[[nm]] <- term_to_pratt[[nm]]
    }
  }
  
  ## ---------------- PI column name ----------------
  if (!is.null(k)) {
    lift_col_name <- "PI_CV"
  } else if (!is.null(pct_train)) {
    lift_col_name <- "PI_Test"
  } else {
    lift_col_name <- "PI"
  }
  
  ## ---------------- helper_table ----------------
  coef_rows <- names(b_hat)
  
  est_vec <- as.numeric(b_hat[coef_rows])
  p_vec   <- as.numeric(p_hat[coef_rows])
  
  cohen_d_vec <- as.numeric(d_hat[coef_rows])
  vif_out_vec <- as.numeric(vif_vec[coef_rows])
  pi_vec      <- as.numeric(pi_by_coef[coef_rows])
  
  helper_table <- data.frame(
    Parameter = coef_rows,
    Coef      = finalfit::round_tidy(est_vec, 2),
    p_value   = paste0(finalfit::round_tidy(p_vec, 3), sapply(p_vec, add_stars)),
    Cohens_d  = finalfit::round_tidy(cohen_d_vec, 2),
    VIF       = finalfit::round_tidy(vif_out_vec, 2),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  helper_table <- cbind(helper_table, PI = finalfit::round_tidy(pi_vec, 3))
  colnames(helper_table)[colnames(helper_table) == "PI"] <- lift_col_name
  
  ## ---------------- summ() + coeftable substitution ----------------
  summ_output <- jtools::summ(stub_lm)
  
  wanted <- c("p_value","Cohens_d","VIF", lift_col_name)
  wanted <- wanted[wanted %in% colnames(helper_table)]
  
  Coef <- align_num_string(helper_table$Coef)
  
  helper_cols <- helper_table[, wanted, drop = FALSE]
  rownames(helper_cols) <- helper_table$Parameter
  
  new_ct <- cbind(Coef = Coef, helper_cols)
  new_ct <- as.data.frame(new_ct, stringsAsFactors = FALSE)
  rownames(new_ct) <- helper_table$Parameter
  
  output_table <- new_ct[, "Coef", drop = FALSE]
  
  if (cohens_d_col && "Cohens_d" %in% colnames(new_ct)) output_table$Cohens_d <- align_num_string(new_ct$Cohens_d)
  if (vif_col     && "VIF"      %in% colnames(new_ct)) output_table$VIF      <- align_num_string(new_ct$VIF)
  if (PI_col      && lift_col_name %in% colnames(new_ct)) output_table[[lift_col_name]] <- align_num_string(new_ct[[lift_col_name]])
  
  output_table$p_value <- new_ct$p_value
  
  # ---- Blank out undefined columns for intercept (match summ_lm behavior) ----
  if ("(Intercept)" %in% rownames(output_table)) {
    int_row <- which(rownames(output_table) == "(Intercept)")
    for (col in intersect(c("VIF", lift_col_name), colnames(output_table))) {
      output_table[int_row, col] <- ""
    }
  }
  
  summ_output[["coeftable"]] <- output_table
  
  rendered <- capture.output(print(summ_output))
  rendered <- gsub("Dependent Variable", "Outcome Variable", rendered, fixed = TRUE)
  
  ## ---------------- insert FE / Clusters in MODEL INFO ----------------
  lab <- .extract_labels(model)
  
  # try Type: first (like before)
  type_idx <- grep("^\\s*Type\\s*:", rendered)
  
  # fallback anchor: MODEL INFO:
  info_idx <- grep("^\\s*MODEL\\s+INFO\\s*:\\s*$", rendered)
  
  anchor <- if (length(type_idx)) type_idx[1] else if (length(info_idx)) info_idx[1] else NA_integer_
  
  if (!is.na(anchor)) {
    indent <- if (anchor < length(rendered)) sub("^([[:space:]]*).*", "\\1", rendered[anchor + 1]) else ""
    ins <- character(0)
    if (!is.null(lab$fixed_effects)) ins <- c(ins, paste0(indent, "Fixed effects: ", lab$fixed_effects))
    if (!is.null(lab$clusters))      ins <- c(ins, paste0(indent, "Clusters: ", lab$clusters))
    if (length(ins)) rendered <- append(rendered, ins, after = anchor)
  }
  
  ## ---------------- HARD REPLACE MODEL FIT BLOCK ----------------
  fit_hdr <- grep("^\\s*MODEL\\s+FIT\\s*:\\s*$", rendered)
  dash_ln <- grep("^\\s*-{5,}\\s*$", rendered)
  
  if (length(fit_hdr) && length(dash_ln)) {
    dash_after <- dash_ln[dash_ln > fit_hdr[1]]
    if (length(dash_after)) {
      cut_at <- dash_after[1]
      indent <- if (fit_hdr[1] < length(rendered)) sub("^([[:space:]]*).*", "\\1", rendered[fit_hdr[1] + 1]) else ""
      
      new_fit <- c(
        paste0(indent, "R² = ", sprintf("%.2f", r2_in)),
        if (is.finite(within_r2)) paste0(indent, "Within R² = ", sprintf("%.2f", within_r2)) else NULL,
        if (!is.finite(within_r2)) paste0(indent, "Adj. R² = ", sprintf("%.2f", adj_in)) else NULL
      )
      
      if (!is.null(k) && is.finite(r2_oos)) {
        lab_cv <- if (isTRUE(model$fes)) "CV Within R²" else "CV R²"
        new_fit <- c(
          new_fit,
          paste0(indent, lab_cv, " = ", sprintf("%.2f", r2_oos)),
          paste0(indent, "Folds = ", k)
        )
      } else if (!is.null(pct_train) && is.finite(r2_oos)) {
        lab_te <- if (isTRUE(model$fes)) "Test Set Within R²" else "Test Set R²"
        new_fit <- c(
          new_fit,
          paste0(indent, lab_te, " = ", sprintf("%.2f", r2_oos)),
          paste0(indent, "Test Set Size = ", sprintf("%.0f%%", (1 - pct_train) * 100))
        )
      }
      
      rendered <- c(
        rendered[1:fit_hdr[1]],
        new_fit,
        "",
        rendered[cut_at:length(rendered)]
      )
    }
  }
  
  ## ---------------- trim noise (same as summ_lm) ----------------
  rendered <- rendered[!grepl("^F\\(", rendered)]
  std_idx  <- grepl("^\\s*Standard\\s+errors\\s*:", rendered, ignore.case = TRUE)
  rendered <- rendered[!std_idx]
  
  ## ---- Force-correct p-value stars post-render (same idea as summ_lm) ----
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
    
    paste0(substr(line, 1, start - 1), replacement, substr(line, start + len, nchar(line)))
  }
  
  rendered <- vapply(rendered, fix_p_line, character(1))
  rendered <- gsub("`+", "", rendered, perl = TRUE)
  
  cat(paste0(rendered, collapse = "\n"), "\n")
  cat("p: *** is <.001, ** is <.01, * is <.05\n",
      "|d|: Trivial <.1 < Tiny <.2 < Small <.5 < Medium <.8 < Large <1.4 < Huge\n",
      sep = "")
  cat("\n")
  
  invisible(NULL)
}