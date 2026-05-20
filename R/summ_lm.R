summ_lm <- function(model, k = NULL, pct_train = NULL,
                    std_est_col = TRUE,
                    sd_est      = FALSE,
                    vif_col      = TRUE,
                    PI_col       = TRUE) {
  
  if (!inherits(model, "lm") && !inherits(model, "lm_robust")) {
    stop("`model` must be a fitted lm or lm_robust object.")
  }
  if (inherits(model, "lm_robust")) return(summ_lm_robust(model, k=k, pct_train=pct_train, sd_est=sd_est, vif_col=vif_col, PI_col=PI_col))
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
  .drop_int <- function(x) x[x != "(Intercept)"]
  
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
  
  # Standardized estimate from coefficient, predictor SD, and outcome SD
  # b_std = b * sd(x) / sd(y)
  std_est_from_fit <- function(coef_names, b, sd_x_vec, sd_y) {
    setNames(b * sd_x_vec / sd_y, coef_names)
  }
  
  # --- aGVIF from glmtoolbox on the ORIGINAL lm at the COEFFICIENT level ---
  .agvif_glmtoolbox_vec <- function(fitted_model) {
    # suppress gvif()'s printed output
    G <- {
      suppressMessages(suppressWarnings(
        capture.output(out <- glmtoolbox::gvif(fitted_model))
      ))
      out
    }
    
    # third column is aGVIF = GVIF^(1/(2*df))
    # square it so it works with standard VIF benchmarks
    a_by_term <- as.numeric(G[, 3])^2
    names(a_by_term) <- rownames(G)
    
    # map term-level aGVIF back to coefficient rows
    X_full   <- stats::model.matrix(fitted_model)
    has_int  <- isTRUE(attr(stats::terms(fitted_model), "intercept") == 1)
    coef_all <- names(stats::coef(fitted_model))
    coef_names <- if (has_int) coef_all[-1] else coef_all
    
    assign_full <- attr(X_full, "assign")
    assign_vec  <- if (has_int) assign_full[-1] else assign_full
    term_labels <- attr(stats::terms(fitted_model), "term.labels")
    coef_terms  <- term_labels[assign_vec]
    
    v <- setNames(a_by_term[coef_terms], coef_names)
    if (has_int) v <- c("(Intercept)" = NA_real_, v)
    v
  }
  
  ## ---------------- evaluated data & names ----------------
  mf        <- model.frame(model)      # evaluated (log(), I(), poly(), etc. already applied)
  tt_eval   <- terms(model)
  y         <- model.response(mf)
  n         <- nrow(mf)
  has_intercept <- isTRUE(attr(tt_eval, "intercept") == 1)
  
  X_full <- model.matrix(tt_eval, mf)
  X      <- if (has_intercept) X_full[, -1, drop = FALSE] else X_full
  
  coef_names  <- colnames(X)
  term_labels <- attr(tt_eval, "term.labels")
  assign_full <- attr(X_full, "assign")
  assign_vec  <- if (has_intercept) assign_full[-1] else assign_full
  coef2term   <- data.frame(variable = coef_names,
                            term     = term_labels[assign_vec],
                            stringsAsFactors = FALSE)
  
  ## ---------------- Standardized estimate from the fitted lm ----------------
  std_est_fun <- function(fitted_lm) {
    sm <- summary(fitted_lm)$coefficients
    rn <- rownames(sm); rn <- rn[rn != "(Intercept)"]
    if (!length(rn)) {
      return(data.frame(variable = character(0),
                        std_est = numeric(0),
                        stringsAsFactors = FALSE))
    }
    mm  <- model.matrix(fitted_lm)
    X_noi <- mm[, rn, drop = FALSE]
    y_fit <- model.response(model.frame(fitted_lm))
    sd_y  <- sd(y_fit)
    sd_x_vec <- apply(X_noi, 2, sd)
    b     <- sm[rn, "Estimate"]
    std_b <- b * sd_x_vec / sd_y
    data.frame(variable = rn, std_est = as.numeric(std_b), stringsAsFactors = FALSE)
  }
  
  std_est_df <- std_est_fun(model)
  
  ## ---------------- VIF (aGVIF) ----------------
  vif_values <- if (ncol(model[["model"]]) > 2) {
    .agvif_glmtoolbox_vec(model)
  } else {
    setNames(rep(NA_real_, length(coef(model))), names(coef(model)))
  }
  
  ## ---------------- model fit R²s ----------------
  r2_in  <- calc_R2(y, fitted(model), has_intercept)
  r2_oos <- NA_real_
  
  ## ----- Train/Test (transform-safe; no lm(fml,...)) -----
  train <- test <- NULL; y_tr <- y_te <- NULL; fit_tr <- NULL
  if (!is.null(pct_train)) {
    if (pct_train <= 0 || pct_train >= 1) stop("`pct_train` must be in (0, 1).")
    idx_tr <- sample(seq_len(n), size = floor(pct_train * n))
    train  <- mf[idx_tr, , drop = FALSE]
    test   <- mf[-idx_tr, , drop = FALSE]
    
    X_tr <- model.matrix(tt_eval, train)
    y_tr <- model.response(train)
    fit_tr <- stats::lm.fit(X_tr, y_tr)
    
    X_te <- model.matrix(tt_eval, test)
    y_te <- model.response(test)
    y_hat_te <- as.vector(X_te %*% fit_tr$coefficients)
    
    r2_oos <- calc_R2(y_te, y_hat_te, has_intercept, mean(y_tr))
  }
  
  ## ----- K-fold CV (transform-safe; no lm(fml,...)) -----
  fold_id <- NULL; y_pred_all <- NULL
  if (!is.null(k)) {
    if (identical(k, "n")) k <- n
    if (k < 2 || k > n) stop("`k` must be an integer between 2 and n.")
    fold_id <- sample(rep(1:k, length.out = n))
    
    y_pred_all <- numeric(n)
    for (fold in seq_len(k)) {
      tr <- mf[fold_id != fold, , drop = FALSE]
      te <- mf[fold_id == fold, , drop = FALSE]
      
      X_tr <- model.matrix(tt_eval, tr)
      y_tr <- model.response(tr)
      fit  <- stats::lm.fit(X_tr, y_tr)
      
      X_te <- model.matrix(tt_eval, te)
      y_pred_all[fold_id == fold] <- as.vector(X_te %*% fit$coefficients)
    }
    r2_oos <- calc_R2(y, y_pred_all, has_intercept)
  }
  
  ## ----- SD of standardized estimate across folds (transform-safe; no lm(fml,...)) -----
  sd_std_e_cv <- rep(NA_real_, length(coef_names))
  names(sd_std_e_cv) <- coef_names
  if (!is.null(k)) {
    std_e_mat <- matrix(NA_real_, nrow = k, ncol = length(coef_names))
    colnames(std_e_mat) <- coef_names
    for (fold in seq_len(k)) {
      tr   <- mf[fold_id != fold, , drop = FALSE]
      X_tr <- model.matrix(tt_eval, tr)
      y_tr <- model.response(tr)
      fit  <- stats::lm.fit(X_tr, y_tr)
      
      sd_y_fold <- sd(y_tr)
      
      all_coef_names <- colnames(X_tr)
      if (has_intercept) {
        keep_idx <- which(all_coef_names != "(Intercept)")
        b_kept   <- fit$coefficients[keep_idx]
        X_kept   <- X_tr[, keep_idx, drop = FALSE]
        names(b_kept) <- all_coef_names[keep_idx]
      } else {
        b_kept   <- fit$coefficients
        X_kept   <- X_tr
        names(b_kept) <- all_coef_names
      }
      
      sd_x_vec   <- apply(X_kept, 2, sd)
      std_e_fold <- std_est_from_fit(names(b_kept), b_kept, sd_x_vec, sd_y_fold)
      
      present <- intersect(coef_names, names(std_e_fold))
      std_e_mat[fold, present] <- std_e_fold[match(present, names(std_e_fold))]
    }
    sd_std_e_cv <- apply(std_e_mat, 2, stats::sd, na.rm = TRUE)
  }
  
  ## ---------------- Pratt-style importance (std_est * correlation) ----------------
  # choose target R² (in-sample vs validation)
  if (!is.null(k) || !is.null(pct_train)) {
    R2_target <- if (!is.na(r2_oos)) r2_oos else r2_in
  } else {
    R2_target <- r2_in
  }
  
  # coefficient-level standardized beta
  mm2      <- stats::model.matrix(model)
  coef_all <- colnames(mm2)
  coef_noi <- coef_all[coef_all != "(Intercept)"]
  
  # raw correlation between each model-matrix column and y
  X_noi <- mm2[, coef_noi, drop = FALSE]
  r_vec <- vapply(seq_along(coef_noi), function(j) {
    xj <- X_noi[, j]
    if (all(is.na(xj)) || stats::sd(xj) == 0) return(0)
    suppressWarnings(stats::cor(xj, y))
  }, numeric(1))
  names(r_vec) <- coef_noi
  
  std_e_vec <- std_est_df$std_est[match(coef_noi, std_est_df$variable)]
  std_e_vec[is.na(std_e_vec)] <- 0
  
  # Pratt: std_est * r. These sum to R²_in by construction in OLS.
  raw_contrib <- std_e_vec * r_vec
  
  pratt_coef_df <- tibble::tibble(
    variable = coef_noi,
    raw      = raw_contrib
  )
  
  # add term labels
  mm_assign <- attr(mm2, "assign")
  mm_terms  <- attr(stats::terms(model), "term.labels")
  term_idx  <- mm_assign[coef_all != "(Intercept)"]
  pratt_coef_df$term <- mm_terms[term_idx]
  
  term_sum <- pratt_coef_df |>
    dplyr::group_by(term) |>
    dplyr::summarise(S = sum(raw, na.rm = TRUE), .groups = "drop")
  
  S_tot <- sum(term_sum$S, na.rm = TRUE)
  
  # Rescale to target R² (matters when validation R² differs from in-sample R²;
  # is essentially a no-op when R2_target == r2_in since S_tot == r2_in by construction).
  if (is.finite(S_tot) && S_tot != 0 && is.finite(R2_target) && !is.na(R2_target)) {
    term_sum <- term_sum |>
      dplyr::mutate(Pratt_term = S / S_tot * R2_target)
  } else {
    term_sum$Pratt_term <- NA_real_
  }
  
  # map term-level Pratt back to coefficient-level (for non-intercept coefs)
  term_to_pratt <- setNames(term_sum$Pratt_term, term_sum$term)
  
  pratt_coef_df$Pratt_coef <- term_to_pratt[pratt_coef_df$term]
  
  # now a coefficient-level vector aligned to coef names (excluding intercept)
  pratt_vec <- setNames(rep(NA_real_, length(coef_all)), coef_all)
  pratt_vec[coef_noi] <- pratt_coef_df$Pratt_coef[match(coef_noi, pratt_coef_df$variable)]
  
  # column name depending on validation mode
  if (!is.null(k)) {
    lift_col_name <- "PI_CV"
  } else if (!is.null(pct_train)) {
    lift_col_name <- "PI_Test"
  } else {
    lift_col_name <- "PI"
  }
  
  ## ---------------- build display table ----------------
  sum_model  <- summary(model)
  coef_rows  <- .drop_int(rownames(sum_model[["coefficients"]]))
  ci_mat     <- confint(model)
  
  mm2       <- stats::model.matrix(model)
  assignV   <- attr(mm2, "assign")
  termLbl   <- attr(stats::terms(model), "term.labels")
  coef_all2  <- names(stats::coef(model))
  assign_no_int <- assignV[coef_all2 != "(Intercept)"]
  map_df <- data.frame(variable = .drop_int(coef_all2),
                       term     = termLbl[assign_no_int],
                       stringsAsFactors = FALSE)
  
  eff_size_df <- dplyr::left_join(std_est_df, map_df, by = "variable")
  
  # attach Pratt at coefficient level
  pratt_df_for_join <- tibble::tibble(
    variable = names(pratt_vec)[names(pratt_vec) != "(Intercept)"],
    !!lift_col_name := as.numeric(pratt_vec[names(pratt_vec) != "(Intercept)"])
  )
  eff_size_df <- dplyr::left_join(eff_size_df, pratt_df_for_join, by = "variable")
  
  est_vec <- sum_model[["coefficients"]][coef_rows, "Estimate"]
  p_vec   <- sum_model[["coefficients"]][coef_rows, "Pr(>|t|)"]
  ci_fmt  <- apply(ci_mat[coef_rows, , drop = FALSE], 1,
                   function(x) sprintf("[%.2f, %.2f]", round(x[1],2), round(x[2],2)))
  idx_eff <- match(coef_rows, eff_size_df$variable)
  std_est_disp_vec <- eff_size_df$std_est[idx_eff]
  vif_vec     <- as.numeric(vif_values[coef_rows])
  sd_std_e_vec <- sd_std_e_cv[match(coef_rows, names(sd_std_e_cv))]
  lift_vec    <- eff_size_df[[lift_col_name]][idx_eff]
  
  base_cols <- data.frame(
    Parameter = coef_rows,
    Coef      = finalfit::round_tidy(est_vec, 2),
    Conf_Int  = ci_fmt,
    p_value   = paste0(finalfit::round_tidy(p_vec, 3), sapply(p_vec, add_stars)),
    Std.Est   = finalfit::round_tidy(std_est_disp_vec, 2),
    SD_Est    = finalfit::round_tidy(sd_std_e_vec, 2),
    VIF       = finalfit::round_tidy(vif_vec, 2),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  helper_table <- cbind(base_cols,
                        PI = finalfit::round_tidy(as.numeric(lift_vec), 3))
  
  # optional intercept row
  if ("(Intercept)" %in% rownames(sum_model[["coefficients"]])) {
    ip <- sum_model[["coefficients"]]["(Intercept)", "Pr(>|t|)"]
    intercept_row <- data.frame(
      Parameter  = "(Intercept)",
      Coef       = finalfit::round_tidy(sum_model[["coefficients"]]["(Intercept)","Estimate"], 2),
      Conf_Int   = sprintf("[%.2f, %.2f]",
                           round(ci_mat["(Intercept)",1],2),
                           round(ci_mat["(Intercept)",2],2)),
      p_value    = paste0(finalfit::round_tidy(ip, 3), add_stars(ip)),
      Std.Est    = "", SD_Est = "", VIF = "", PI = "",
      stringsAsFactors = FALSE, check.names = FALSE
    )
    helper_table <- rbind(intercept_row, helper_table)
  }
  rownames(helper_table) <- NULL
  
  dep_var <- deparse(formula(model)[[2]])
  internal_pi_col <- "PI"
  nice_map <- c(
    Parameter = paste0("DV: ", dep_var),
    Coef      = "Coef",
    Conf_Int  = "CI_95%",
    p_value   = "p_value",
    Std.Est   = "Std.Est",
    SD_Est    = "SD(Est.)",
    VIF       = "VIF"
  )
  colnames_df <- vapply(names(helper_table), function(n) {
    if (n == internal_pi_col) return(lift_col_name)
    if (n %in% names(nice_map)) return(nice_map[[n]])
    n
  }, character(1), USE.NAMES = FALSE)
  colnames(helper_table) <- colnames_df
  
  ## ---------------- jtools::summ and print tweaks ----------------
  summ_output <- jtools::summ(model)
  
  old_rn <- rownames(summ_output[["coeftable"]])
  Coef_raw <- summ_output[["coeftable"]][, 1]
  Coef     <- align_num_string(Coef_raw)
  
  wanted <- c("p_value","Std.Est","SD(Est.)","VIF", lift_col_name)
  wanted <- wanted[wanted %in% colnames(helper_table)]
  
  helper_rows <- match(c("(Intercept)", old_rn), helper_table[[1]])
  helper_cols <- helper_table[ helper_rows, wanted, drop = FALSE ]
  rownames(helper_cols) <- make.unique(c("(Intercept)", old_rn))
  
  new_ct <- cbind(Coef = Coef, helper_cols[ make.unique(old_rn), , drop = FALSE ])
  new_ct <- as.data.frame(new_ct, stringsAsFactors = FALSE)
  rownames(new_ct) <- make.unique(old_rn)
  
  output_table <- new_ct[, 1, drop = FALSE]
  if (std_est_col && "Std.Est" %in% colnames(new_ct)) output_table$`Std.Est` <- new_ct$`Std.Est`
  if (std_est_col && sd_est && !is.null(k) && "SD(Est.)" %in% colnames(new_ct)) output_table$`SD_Est` <- new_ct$`SD_Est`
  if (vif_col && "VIF" %in% colnames(new_ct)) output_table$`VIF` <- new_ct$`VIF`
  if (PI_col && lift_col_name %in% colnames(new_ct)) {
    output_table[[lift_col_name]] <- align_num_string(new_ct[[lift_col_name]])
  }
  if ("Std.Est" %in% colnames(output_table)) {
    output_table$Std.Est <- align_num_string(output_table$Std.Est)
  }
  if ("SD(Est.)" %in% colnames(output_table)) {
    output_table$SD_Est <- align_num_string(output_table$SD_Est)
  }
  output_table$`p_value` <- new_ct$`p_value`
  summ_output[["coeftable"]] <- output_table
  
  rendered <- capture.output(print(summ_output))
  rendered <- gsub("Dependent Variable", "Outcome Variable", rendered, fixed = TRUE)
  
  label <- NULL; val <- NULL; extra <- NULL
  if (!is.null(k)) {
    label <- "CV R²";      val <- sprintf("%.2f", r2_oos)
    extra <- function(indent) paste0(indent, "Folds = ", k)
  } else if (!is.null(pct_train)) {
    label <- "Test Set R²"; val <- sprintf("%.2f", r2_oos)
    test_pct <- (1 - pct_train) * 100
    extra <- function(indent) paste0(indent, "Test Set Size = ", sprintf("%.0f%%", test_pct))
  }
  if (!is.null(label)) {
    adj_idx <- grep("^\\s*Adj", rendered, ignore.case = TRUE)
    if (length(adj_idx)) {
      indent <- sub("^([[:space:]]*).*", "\\1", rendered[adj_idx[1]])
      rendered[adj_idx[1]] <- paste0(indent, label, " = ", val)
      rendered <- append(rendered, extra(indent), after = adj_idx[1])
    }
  }
  
  # trim noise
  rendered <- rendered[!grepl("^F\\(", rendered)]
  std_idx  <- grepl("^\\s*Standard\\s+errors\\s*:", rendered, ignore.case = TRUE)
  rendered <- rendered[!std_idx]
  
  ## ---- Force-correct p-value stars post-render (p always "0.xxx") ----
  p_stars <- function(pv_num) {
    if (is.na(pv_num)) return(" ")
    if (pv_num < 0.001) " ***"
    else if (pv_num < 0.01) " **"
    else if (pv_num < 0.05) " *"
    else " "
  }
  
  fix_p_line <- function(line) {
    # find rightmost token "0.xxx"
    m <- gregexpr("\\b0\\.\\d{3}\\b(\\s\\*{0,3})?", line, perl = TRUE)
    pos <- m[[1]]
    if (length(pos) == 1 && pos[1] == -1) return(line)
    
    i <- tail(which(pos > 0), 1)
    start <- pos[i]
    len   <- attr(m[[1]], "match.length")[i]
    
    tok <- substr(line, start, start + len - 1)
    # strip any existing space+stars after the number
    num_str <- sub("\\s\\*{0,3}$", "", tok)
    pv <- suppressWarnings(as.numeric(num_str))
    
    stars <- p_stars(pv)
    # rebuild: number + correct stars
    replacement <- paste0(num_str, stars)
    
    paste0(substr(line, 1, start - 1), replacement, substr(line, start + len, nchar(line)))
  }
  
  rendered <- vapply(rendered, fix_p_line, character(1))
  rendered <- gsub("`+", "", rendered, perl = TRUE)
  
  cat(paste0(rendered, collapse = "\n"), "\n")
  cat("p: *** is <.001, ** is <.01, * is <.05\n",
      "|Std.Est|: Trivial <.1 < Small <.3 < Medium <.5 < Large\n",
      sep = "")
  cat("\n")
}
