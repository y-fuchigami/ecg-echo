proportion_ci <- function(x, n, method = "wilson") {
  if (is.na(n) || n == 0) {
    return(c(NA_real_, NA_real_))
  }
  if (method == "exact") {
    interval <- binom.test(x, n)$conf.int
    return(unname(interval))
  }
  ci <- binom::binom.confint(x, n, methods = "wilson")
  c(ci$lower, ci$upper)
}

log_ci_ratio <- function(est, se_log) {
  log_est <- log(est)
  margin <- c(-1, 1) * 1.96 * se_log
  exp(log_est + margin)
}

diagnostic_one <- function(data, exposure, outcome,
                           label_exposure = exposure, label_outcome = outcome,
                           ci_method = "wilson") {
  df <- data
  df$exp_num <- to01(df[[exposure]])
  df$out_num <- to01(df[[outcome]])
  keep <- !is.na(df$exp_num) & !is.na(df$out_num)
  df <- df[keep]

  tp <- sum(df$exp_num == 1 & df$out_num == 1)
  fp <- sum(df$exp_num == 1 & df$out_num == 0)
  fn <- sum(df$exp_num == 0 & df$out_num == 1)
  tn <- sum(df$exp_num == 0 & df$out_num == 0)

  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)

  sens_ci <- proportion_ci(tp, tp + fn, ci_method)
  spec_ci <- proportion_ci(tn, tn + fp, ci_method)
  ppv_ci <- proportion_ci(tp, tp + fp, ci_method)
  npv_ci <- proportion_ci(tn, tn + fn, ci_method)

  lr_pos <- sens / (1 - spec)
  lr_neg <- (1 - sens) / spec

  if (tp > 0 && fn > 0 && fp > 0 && tn > 0) {
    se_pos <- sqrt(1 / tp - 1 / (tp + fn) + 1 / fp - 1 / (fp + tn))
    se_neg <- sqrt(1 / fn - 1 / (tp + fn) + 1 / tn - 1 / (fp + tn))
    lr_pos_ci <- log_ci_ratio(lr_pos, se_pos)
    lr_neg_ci <- log_ci_ratio(lr_neg, se_neg)
  } else {
    lr_pos_ci <- c(NA_real_, NA_real_)
    lr_neg_ci <- c(NA_real_, NA_real_)
  }

  data.table::data.table(
    exposure = label_exposure, outcome = label_outcome,
    n = nrow(df), tp = tp, fp = fp, fn = fn, tn = tn,
    sensitivity = sens, sensitivity_low = sens_ci[1], sensitivity_high = sens_ci[2],
    specificity = spec, specificity_low = spec_ci[1], specificity_high = spec_ci[2],
    ppv = ppv, ppv_low = ppv_ci[1], ppv_high = ppv_ci[2],
    npv = npv, npv_low = npv_ci[1], npv_high = npv_ci[2],
    lr_pos = lr_pos, lr_pos_low = lr_pos_ci[1], lr_pos_high = lr_pos_ci[2],
    lr_neg = lr_neg, lr_neg_low = lr_neg_ci[1], lr_neg_high = lr_neg_ci[2]
  )
}
