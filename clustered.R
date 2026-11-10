fit_clustered_logit <- function(data, outcome) {
  col <- paste0(outcome, "_num")
  keep <- !is.na(data[[col]])
  keep <- keep & !is.na(data$rad_or_rvh_num)
  keep <- keep & !is.na(data$subject_id)
  df <- data[keep]

  covars <- paste(covars_full_model, collapse = " + ")
  form <- paste(col, "~ rad_or_rvh_num +", covars)
  form <- as.formula(form)
  fit <- glm(form, data = df, family = binomial())

  vcov <- sandwich::vcovCL(fit, cluster = df$subject_id, type = "HC0")
  coefs <- lmtest::coeftest(fit, vcov. = vcov)
  estimate <- coefs["rad_or_rvh_num", "Estimate"]
  se <- coefs["rad_or_rvh_num", "Std. Error"]

  used <- rep(TRUE, nrow(df))
  used[fit$na.action] <- FALSE
  df <- df[used]

  data.table::data.table(
    outcome = outcome,
    n_pairs = nrow(df),
    n_patients = data.table::uniqueN(df$subject_id),
    n_events = sum(df[[col]] == 1),
    n_exposed = sum(df$rad_or_rvh_num == 1),
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * se),
    CI_high = exp(estimate + 1.96 * se),
    p_value = coefs["rad_or_rvh_num", "Pr(>|z|)"]
  )
}
