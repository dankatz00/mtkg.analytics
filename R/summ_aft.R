summ_aft <- function(fit,
                     time_horizon,
                     center = FALSE,
                     k = NULL,
                     pct_train = NULL,
                     PI_col = TRUE,
                     pct_change_col = TRUE,
                     vif_col = TRUE) {
  
  stopifnot(inherits(fit, "survreg"))
  if (is.null(time_horizon) || !is.finite(time_horizon) || time_horizon <= 0) {
    stop("`time_horizon` must be a positive number on the same time scale as Surv(time, status).")
  }
  set.seed(123)
  
  # -- silent, robust package ensure ------------------------------
  pkgs <- c("survival","jtools","dplyr","tibble","finalfit","glmtoolbox")
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
  
  # ---------------- helpers (match summ_logit style) ----------------
  tjur_r2 <- function(y, p) { y <- as.numeric(y); mean(p[y==1]) - mean(p[y==0]) }
  
  align_num_string <- function(x, digits = 2) {    # fixed-width numeric
    num <- suppressWarnings(as.numeric(x))
    s   <- ifelse(is.na(num), "", sprintf(paste0("%.", digits, "f"), num))
    s   <- ifelse(startsWith(s,"-"), s, paste0("\u2007", s))   # figure-space for positives
    w   <- max(nchar(s, type="width"), na.rm = TRUE)
    paste0(strrep("\u2007", pmax(0L, w - nchar(s, type="width"))), s)
  }
  
  fw <- function(chars) {                           # fixed-width strings (e.g., p-values)
    w <- max(nchar(chars, type="width"), na.rm = TRUE)
    paste0(strrep("\u2007", w - nchar(chars, type="width")), chars)
  }
  
  add_stars_placeholder <- function(p) {            # placeholders so the dashed rule length doesn't shrink
    if (is.na(p)) return("\u2007")
    if (p < 0.001) return("\u2007````")
    if (p < 0.01)  return("\u2007```")
    if (p < 0.05)  return("\u2007``")
    "\u2007"
  }
  
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
    
    replacement <- paste0(num_str, p_stars(pv))
    paste0(substr(line, 1, start - 1), replacement,
           substr(line, start + len, nchar(line)))
  }
  
  # Extract the actual time/status variable names from Surv(time,status)
  .surv_lhs_names <- function(fml) {
    lhs <- fml[[2]]
    if (!is.call(lhs) || as.character(lhs[[1]]) != "Surv") {
      stop("Model formula must have Surv(time, status) on the left-hand side.")
    }
    if (length(lhs) < 3) stop("Surv() must have at least time and status.")
    tnm <- as.character(lhs[[2]])
    snm <- as.character(lhs[[3]])
    list(time = tnm, status = snm)
  }
  
  # --- aGVIF from glmtoolbox on the OLS stub, returned at the COEFFICIENT level ---
  .agvif_glmtoolbox_vec <- function(lm_stub) {
    suppressMessages(suppressWarnings(capture.output(G <- glmtoolbox::gvif(lm_stub))))
    
    # third column is aGVIF = GVIF^(1/(2*df)) squared (so same scale as regular VIF)
    a_by_term <- as.numeric(G[, 3])^2
    names(a_by_term) <- rownames(G)
    
    X_full  <- stats::model.matrix(lm_stub)
    has_int <- isTRUE(attr(stats::terms(lm_stub), "intercept") == 1)
    
    coef_all   <- names(stats::coef(lm_stub))
    coef_names <- if (has_int) coef_all[-1] else coef_all
    
    assign_full <- attr(X_full, "assign")
    assign_vec  <- if (has_int) assign_full[-1] else assign_full
    term_labels <- attr(stats::terms(lm_stub), "term.labels")
    coef_terms  <- term_labels[assign_vec]
    
    v <- setNames(a_by_term[coef_terms], coef_names)
    if (has_int) v <- c("(Intercept)" = NA_real_, v)
    v
  }
  
  # P(fail by t0 | X) from survreg using psurvreg on the log-time scale
  pred_fail_prob <- function(fitted_survreg, newdata, t0) {
    lp <- stats::predict(fitted_survreg, newdata = newdata, type = "lp")
    survival::psurvreg(q = log(t0),
                       mean = lp,
                       scale = fitted_survreg$scale,
                       distribution = fitted_survreg$dist)
  }
  
  build_y_t0 <- function(time, status, t0) {
    y <- rep(NA_integer_, length(time))
    known <- rep(FALSE, length(time))
    
    idx1 <- (time <= t0) & (status == 1)  # failed by t0
    y[idx1] <- 1L
    known[idx1] <- TRUE
    
    idx0 <- (time > t0)                   # known survived past t0
    y[idx0] <- 0L
    known[idx0] <- TRUE
    
    list(y = y, known = known)
  }
  
  # --------- VALIDATION ----------
  validate_tjur_aft <- function(f, df_refit, y_sc, t0, k, pct_train) {
    n_sc <- nrow(df_refit)
    if (n_sc < 5) stop("Too few known-by-horizon observations at this time_horizon.")
    
    tjur_fun <- function(y_vec, p_vec) {
      y_vec <- as.numeric(y_vec)
      mean(p_vec[y_vec == 1], na.rm = TRUE) - mean(p_vec[y_vec == 0], na.rm = TRUE)
    }
    
    if (!is.null(k)) {
      if (is.character(k)) {
        if (identical(k, "n")) k <- n_sc else stop("`k` must be an integer or 'n'.")
      }
      k <- as.integer(k)
      if (k < 2L || k > n_sc) stop("`k` must be between 2 and n (or 'n') after filtering by horizon.")
      
      folds <- sample(rep_len(seq_len(k), n_sc))
      p_all <- numeric(n_sc)
      
      for (fold in seq_len(k)) {
        tr <- df_refit[folds != fold, , drop = FALSE]
        te <- df_refit[folds == fold, , drop = FALSE]
        
        mi <- survival::survreg(stats::formula(f), data = tr, dist = f$dist)
        p_all[folds == fold] <- pred_fail_prob(mi, te, t0)
      }
      
      val <- tjur_fun(y_sc, p_all)
      return(list(kind = "CV", val = val, k = k, n_test = NA_integer_, n_scored = n_sc))
    }
    
    if (!is.null(pct_train) && pct_train > 0 && pct_train < 1) {
      ntr <- floor(pct_train * n_sc)
      idx <- sample.int(n_sc, ntr)
      tr  <- df_refit[idx, , drop = FALSE]
      te  <- df_refit[-idx, , drop = FALSE]
      
      mi   <- survival::survreg(stats::formula(f), data = tr, dist = f$dist)
      p_te <- pred_fail_prob(mi, te, t0)
      y_te <- y_sc[-idx]
      
      val <- tjur_fun(y_te, p_te)
      return(list(kind = "Holdout", val = val, k = NA_integer_, n_test = nrow(te), n_scored = n_sc))
    }
    
    NULL
  }
  
  # ---------------- evaluated data ----------------
  mf <- stats::model.frame(fit)
  
  if (center) {
    
    # response column name (exclude from centering)
    resp_idx  <- attr(stats::terms(fit), "response")
    resp_name <- names(mf)[resp_idx]
    
    is_centerable <- vapply(names(mf), function(nm) {
      x <- mf[[nm]]
      is.numeric(x) &&
        nm != resp_name &&
        length(unique(x[!is.na(x)])) > 2
    }, logical(1))
    
    mf[is_centerable] <- lapply(mf[is_centerable], function(x) x - mean(x, na.rm = TRUE))
    
    # REFIT with a fresh formula environment
    f <- stats::formula(fit)
    environment(f) <- environment()
    
    fit <- survival::survreg(f, data = mf, dist = fit$dist)
    
    # refresh mf for downstream alignment
    mf <- stats::model.frame(fit)
  }
  
  yS <- stats::model.response(mf)
  if (!inherits(yS, "Surv")) stop("Response must be Surv(time, status).")
  time   <- yS[, "time"]
  status <- yS[, "status"]
  
  n_total <- nrow(mf)
  
  yt <- build_y_t0(time, status, time_horizon)
  known <- yt$known
  y_all <- yt$y
  pct_unknown <- mean(!known) * 100  # kept internal (still useful for debugging)
  
  mf_sc <- mf[known, , drop = FALSE]
  y_sc  <- y_all[known]
  if (length(y_sc) < 5) stop("Too few known-by-horizon observations at this `time_horizon`.")
  
  # Build a refit-ready scored data frame that explicitly includes time/status with correct names
  nm <- .surv_lhs_names(stats::formula(fit))
  df_refit_sc <- mf_sc
  df_refit_sc[[nm$time]]   <- time[known]
  df_refit_sc[[nm$status]] <- status[known]
  
  # ---------------- in-sample Tjur at t0 ----------------
  p_hat <- pred_fail_prob(fit, df_refit_sc, time_horizon)
  tjur_train <- tjur_r2(y_sc, p_hat)
  
  # adjusted Tjur (OLS-style adjustment, using scored N)
  p_num <- length(stats::coef(fit)) - 1L
  n_sc  <- length(y_sc)
  tjur_adj <- 1 - ((1 - tjur_train) * (n_sc - 1) / (n_sc - p_num - 1))
  
  # ---------------- optional validation (CV/test) ----------------
  vinfo <- validate_tjur_aft(fit, df_refit_sc, y_sc, time_horizon, k, pct_train)
  
  if (is.null(vinfo)) {
    R2_target <- tjur_train
    lift_col_name <- "PI"
  } else if (identical(vinfo$kind, "Holdout")) {
    R2_target <- vinfo$val
    lift_col_name <- "PI_Test"
  } else {
    R2_target <- vinfo$val
    lift_col_name <- "PI_CV"
  }
  
  # ---------------- coeff table from survreg ----------------
  sm <- as.data.frame(summary(fit)$table)
  sm$Parameter <- rownames(sm)
  rownames(sm) <- NULL
  
  # drop Log(scale) row from display/PI
  sm <- sm[sm$Parameter != "Log(scale)", , drop = FALSE]
  
  beta_raw <- sm$Value
  p_raw    <- sm$p
  
  exp_beta <- exp(beta_raw)
  pct_chg  <- (exp_beta - 1) * 100
  d_raw    <- beta_raw / fit$scale
  
  # blank intercept effect columns
  is_int <- sm$Parameter == "(Intercept)"
  pct_chg[is_int] <- NA_real_
  d_raw[is_int]   <- NA_real_
  
  # ---------------- Pratt-style PI at horizon (term-level scaling, coef-level display) ----------------
  X_mm <- stats::model.matrix(stats::terms(fit), df_refit_sc)
  coef_all <- colnames(X_mm)
  coef_noi <- coef_all[coef_all != "(Intercept)"]
  
  X_noi <- X_mm[, coef_noi, drop = FALSE]
  
  r_vec <- vapply(seq_along(coef_noi), function(j) {
    xj <- X_noi[, j]
    if (all(is.na(xj)) || stats::sd(xj) == 0) return(0)
    suppressWarnings(stats::cor(xj, y_sc))
  }, numeric(1))
  names(r_vec) <- coef_noi
  
  d_named <- d_raw
  names(d_named) <- sm$Parameter
  
  d_vec <- d_named[coef_noi]
  d_vec[is.na(d_vec)] <- 0
  
  raw_contrib <- d_vec * r_vec
  
  mm_assign   <- attr(X_mm, "assign")
  term_labels <- attr(stats::terms(fit), "term.labels")
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
  
  # ---------------- build a jtools skeleton via OLS stub (printing only) ----------------
  rhs <- stats::formula(fit)[[3]]
  rhs_txt <- paste(deparse(rhs), collapse = "")
  lm_formula <- stats::as.formula(paste0("log(time) ~ ", rhs_txt))
  
  df_stub <- mf
  df_stub$time <- time
  df_stub$status <- status
  
  df_events <- df_stub[df_stub$status == 1, , drop = FALSE]
  if (nrow(df_events) < 3) stop("Need at least 3 events (status==1) to build the print skeleton.")
  
  lm_stub <- stats::lm(lm_formula, data = df_events)
  s_obj   <- jtools::summ(lm_stub, digits = 3)
  
  # ---------------- aGVIF / VIF from the OLS stub ----------------
  vifs <- NULL
  if (isTRUE(vif_col)) {
    if (ncol(lm_stub[["model"]]) > 2) {
      vifs <- .agvif_glmtoolbox_vec(lm_stub)  # coefficient-level, includes (Intercept)=NA
    } else {
      vifs <- setNames(rep(NA_real_, length(stats::coef(lm_stub))), names(stats::coef(lm_stub)))
    }
  }
  
  # ---------------- build coeftable (row-aligned to skeleton) ----------------
  rn <- rownames(s_obj[["coeftable"]])
  idx <- match(rn, sm$Parameter)
  if (anyNA(idx)) {
    bad <- rn[is.na(idx)]
    stop("Could not align coefficient rows to survreg coefficients: ", paste(bad, collapse = ", "))
  }
  
  exp_col <- exp_beta[idx]
  pc_col  <- pct_chg[idx]
  d_col   <- d_raw[idx]
  p_col   <- p_raw[idx]
  
  ct <- data.frame(matrix("", nrow = length(rn), ncol = 0), check.names = FALSE)
  rownames(ct) <- rn
  
  ct[["exp(Beta)"]] <- align_num_string(exp_col, 2)
  if (isTRUE(pct_change_col)) ct[["% Change"]] <- align_num_string(pc_col, 2)
  ct[["Cohens_d"]] <- align_num_string(d_col, 2)
  
  if (isTRUE(vif_col)) {
    vif_here <- as.numeric(vifs[rn])
    ct[["VIF"]] <- align_num_string(vif_here, 2)
  }
  
  if (isTRUE(PI_col)) {
    PI_vals <- as.numeric(pratt_vec[rn])
    PI_vals[rn == "(Intercept)"] <- NA_real_
    ct[[lift_col_name]] <- align_num_string(PI_vals, 2)
  }
  
  # p-values: fixed 3dp + placeholder stars, then post-fix
  p_str <- ifelse(is.na(p_col), "", sprintf("%.3f", round(p_col, 3)))
  p_str <- ifelse(p_str == "", "", paste0(p_str, vapply(p_col, add_stars_placeholder, character(1))))
  ct[["p_value"]] <- fw(p_str)
  
  # blank undefined columns for intercept
  if ("(Intercept)" %in% rownames(ct)) {
    int_row <- which(rownames(ct) == "(Intercept)")
    for (col in intersect(c("% Change","Cohens_d","VIF", lift_col_name), colnames(ct))) {
      ct[int_row, col] <- ""
    }
  }
  
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
  
  s_obj[["coeftable"]] <- as.matrix(ct)
  
  # ---------------- render + edits ----------------
  L <- capture.output(print(s_obj))
  
  # Outcome variable fixed text (avoid leaking log(time) from stub)
  L <- gsub("^\\s*Outcome Variable:.*$",
            "Outcome Variable: Time-to-Failure",
            L)
  L <- gsub("^\\s*Dependent Variable:.*$",
            "Outcome Variable: Time-to-Failure",
            L)
  
  # Type line
  L <- sub("^Type:.*$",
           paste0("Type: Accelerated Failure Time (", fit$dist, ")"),
           L)
  
  # Observations: total N
  obs_idx <- grep("^\\s*Observations\\s*:\\s*", L)
  if (length(obs_idx)) {
    L[obs_idx[1]] <- sub("^\\s*Observations\\s*:\\s*\\d+",
                         paste0("Observations: ", n_total),
                         L[obs_idx[1]])
  }
  
  # ---- Insert horizon info under MODEL INFO (after Type) ----
  type_idx <- grep("^Type:", L)
  if (length(type_idx) == 1) {
    indent <- if (type_idx < length(L)) sub("^([[:space:]]*).*", "\\1", L[type_idx + 1]) else ""
    ins <- c(
      paste0(indent, "Time Horizon = ", time_horizon)
    )
    L <- append(L, ins, after = type_idx[1])
  }
  
  # Trim noise
  L <- L[!grepl("^F\\(", L)]
  L <- L[!grepl("^Standard errors:", L)]
  
  # Replace MODEL FIT section with Tjur R² (+ validation if requested) -- no t0 or %unknown lines
  fit_i <- which(grepl("^MODEL FIT:", L))
  if (length(fit_i) == 1) {
    after <- fit_i + 1L
    while (after <= length(L) && nzchar(trimws(L[after]))) after <- after + 1L
    
    r2 <- "\u00B2"
    
    if (is.null(vinfo)) {
      fit_lines <- c(
        "MODEL FIT:",
        sprintf("Tjur's R%s = %.2f", r2, tjur_train),
        sprintf("Adjusted Tjur's R%s = %.2f", r2, tjur_adj)
      )
    } else if (identical(vinfo$kind, "Holdout")) {
      fit_lines <- c(
        "MODEL FIT:",
        sprintf("Tjur's R%s = %.2f", r2, tjur_train),
        sprintf("Tjur's R%s (test) = %.2f", r2, vinfo$val),
        sprintf("Test set size = %.1f%%", (vinfo$n_test / vinfo$n_scored) * 100)
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
  
  # Force-correct p-value stars post-render + strip placeholder backticks
  L <- vapply(L, fix_p_line, character(1))
  L <- gsub("`+", "", L, perl = TRUE)
  
  cat(paste0(paste(L, collapse = "\n"), "\n"))
  cat("p: *** is <.001, ** is <.01, * is <.05\n",
      "|d|: Trivial <.1 < Tiny <.2 < Small <.5 < Medium <.8 < Large <1.4 < Huge\n\n", sep = "")
  
  invisible(s_obj)
}