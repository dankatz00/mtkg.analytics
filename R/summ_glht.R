summ_glht <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    
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
    
    # Compute Standardized Estimates
    mf <- model.frame(x$model)
    y <- model.response(mf)
    X <- model.matrix(x$model)
    K <- x$linfct
    
    sd_y <- sd(y)
    sd_cx <- apply(K, 1, function(k) sd(X %*% k))
    pq$Std.Est <- (pq$coefficients / sd_y) * sd_cx
    
    # Building Coef Table
    mtests <- cbind(round(pq$coefficients, 2), round(pq$Std.Est, 2), pq$pvalues)
    error <- attr(pq$pvalues, "error")
    colnames(mtests) <- c("Est.", "Std.Est", "p")
    type <- pq$type
    
    ### print p values according to simulation precision
    if (!is.null(error) && error > .Machine$double.eps) {
        sig <- which.min(abs(1 / error - (10^(1:10))))
        sig <- 1 / (10^sig)
    } else {
        sig <- .Machine$double.eps
    }
    cat("Linear Hypotheses:\n")
    alt <- switch(x$alternative,
                  "two.sided" = "==", "less" = ">=", "greater" = "<=")
    rownames(mtests) <- paste(rownames(mtests), alt, x$rhs)
    printCoefmat(mtests, digits = digits,
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
