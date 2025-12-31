predict_aft_fail_rate <- function(fit, newdata, time_horizon, return_type = c("prob","data.frame")) {
  stopifnot(inherits(fit, "survreg"))
  if (missing(newdata) || is.null(newdata)) stop("`newdata` is required.")
  if (!is.data.frame(newdata)) stop("`newdata` must be a data.frame.")
  if (is.null(time_horizon) || !is.finite(time_horizon) || time_horizon <= 0) {
    stop("`time_horizon` must be a positive number.")
  }
  return_type <- match.arg(return_type)
  
  # Predicted linear predictor = E[log(T)] (location) for newdata
  lp <- stats::predict(fit, newdata = newdata, type = "lp")
  
  # For survreg Weibull/lognormal/loglogistic, distribution CDF on log-time scale:
  # P(T <= t0 | X) = psurvreg(log(t0), mean = lp, scale = fit$scale, distribution = fit$dist)
  p_fail <- survival::psurvreg(
    q = log(time_horizon),
    mean = lp,
    scale = fit$scale,
    distribution = fit$dist
  )
  
  if (return_type == "prob") return(as.numeric(p_fail))
  
  out <- newdata
  out$pred_fail_rate <- as.numeric(p_fail)
  out
}