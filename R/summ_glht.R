summ_glht <- function(x, digits = 2, ...) {
    
    x <- summary(x, ...)
    
    cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
    if (!is.null(x$type))
        cat("Multiple Comparisons of Means:", x$type, "Contrasts\n\n\n")
    call <- if (isS4(x$model)) x$model@call else x$model$call
    if (!is.null(call)) {
        cat("Fit: ")
        print(call)
        cat("\n")
    }
    
    pq <- x$test
    
    # Compute Standardized Estimates as deviation from null in SD units
    mf <- model.frame(x$model)
    y  <- model.response(mf)
    X  <- model.matrix(x$model)
    K  <- x$linfct
    sd_y <- sd(y)
    rhs  <- x$rhs
    
    # Identify pure main-effect dummy columns of X
    tt <- terms(x$model)
    assign_vec  <- attr(X, "assign")
    term_labels <- attr(tt, "term.labels")
    has_int     <- attr(tt, "intercept") == 1
    
    is_pure_dummy_col <- vapply(seq_len(ncol(X)), function(j) {
        if (has_int && j == 1) return(FALSE)
        col_term <- term_labels[assign_vec[j]]
        if (is.na(col_term)) return(FALSE)
        if (grepl(":", col_term)) return(FALSE)
        all(X[, j] %in% c(0, 1))
    }, logical(1))
    
    deviation <- pq$coefficients - rhs
    
    std_est <- numeric(nrow(K))
    for (i in seq_len(nrow(K))) {
        k_row <- K[i, ]
        touched <- which(k_row != 0)
        touched_no_int <- if (has_int) setdiff(touched, 1) else touched
        
        is_intercept_only <- has_int && length(touched_no_int) == 0 && 1 %in% touched
        is_categorical    <- length(touched_no_int) > 0 && all(is_pure_dummy_col[touched_no_int])
        
        if (is_intercept_only || is_categorical) {
            std_est[i] <- deviation[i] / sd_y
        } else {
            sd_cx <- sd(X %*% k_row)
            std_est[i] <- (deviation[i] / sd_y) * sd_cx
        }
    }
    
    mtests <- cbind(round(deviation, digits),
                    round(std_est, digits),
                    round(pq$pvalues, 3))
    error <- attr(pq$pvalues, "error")
    colnames(mtests) <- c("Est.", "Std.Est", "p")
    type <- pq$type
    
    if (!is.null(error) && error > .Machine$double.eps) {
        sim_sig <- which.min(abs(1 / error - (10^(1:10))))
        sim_sig <- 1 / (10^sim_sig)
        sig <- max(sim_sig, 1e-3)
    } else {
        sig <- 1e-3
    }
    
    cat("Linear Hypotheses:\n")
    alt <- switch(x$alternative,
                  "two.sided" = "==", "less" = ">=", "greater" = "<=")
    
    # Build row labels in deviation-from-zero form: "coef +/- rhs alt 0"
    build_label <- function(coef_label, rhs_val) {
        if (rhs_val == 0) {
            paste(coef_label, alt, 0)
        } else {
            sign_str <- if (rhs_val > 0) "-" else "+"
            paste0(coef_label, " ", sign_str, " ", abs(rhs_val), " ", alt, " 0")
        }
    }
    rownames(mtests) <- mapply(build_label,
                               rownames(mtests), rhs,
                               USE.NAMES = FALSE)
    
    printCoefmat(mtests, digits = digits + 3,
                 has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
    
    switch(type,
        "univariate" = cat("(Univariate p values reported)"),
        "single-step" = cat("(Adjusted p values reported -- single-step method)"),
        "Shaffer" = cat("(Adjusted p values reported -- Shaffer method)"),
        "Westfall" = cat("(Adjusted p values reported -- Westfall method)"),
        cat("(Adjusted p values reported --", type, "method)")
    )
    cat("\n\n")
    invisible(x)
}
