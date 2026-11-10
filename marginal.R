mean_prob <- function(fit, df, value) {
  df$rad_or_rvh_num <- rep(value, nrow(df))
  pred <- predict(fit, newdata = df, type = "response")
  mean(pred, na.rm = TRUE)
}

fit_marginal_probs <- function(data, outcome, B = 300) {
  col <- paste0(outcome, "_num")
  keep <- !is.na(data[[col]])
  keep <- keep & !is.na(data$rad_or_rvh_num)
  keep <- keep & !is.na(data$subject_id)
  df <- data[keep]

  covars <- paste(covars_full_model, collapse = " + ")
  form <- paste(col, "~ rad_or_rvh_num +", covars)
  form <- as.formula(form)
  fit <- glm(form, data = df, family = binomial())

  used <- rep(TRUE, nrow(df))
  used[fit$na.action] <- FALSE
  model_df <- df[used]
  p0 <- mean_prob(fit, df, 0L)
  p1 <- mean_prob(fit, df, 1L)

  set.seed(123)
  ids <- as.character(unique(df$subject_id))
  rows <- split(seq_len(nrow(df)), as.character(df$subject_id))
  boot <- data.table::data.table(
    p0 = numeric(), p1 = numeric(), rd = numeric(), rr = numeric()
  )

  for (i in seq_len(B)) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    sampled_rows <- rows[sampled]
    index <- unlist(sampled_rows, use.names = FALSE)
    boot_df <- df[index]
    if (length(unique(boot_df[[col]])) < 2 ||
        length(unique(boot_df$rad_or_rvh_num)) < 2) {
      next
    }

    boot_fit <- tryCatch(
      glm(form, data = boot_df, family = binomial()),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) {
      next
    }

    p0b <- mean_prob(boot_fit, boot_df, 0L)
    p1b <- mean_prob(boot_fit, boot_df, 1L)
    boot <- rbind(boot, data.table::data.table(
      p0 = p0b, p1 = p1b, rd = p1b - p0b, rr = p1b / p0b
    ))
  }

  data.table::data.table(
    outcome = outcome,
    outcome_label = outcome_labels[[outcome]],
    n_pairs = nrow(model_df),
    n_patients = data.table::uniqueN(model_df$subject_id),
    n_events = sum(model_df[[col]] == 1),
    adjusted_probability_unexposed = p0,
    adjusted_probability_exposed = p1,
    adjusted_risk_difference = p1 - p0,
    adjusted_risk_ratio = p1 / p0,
    p0_low = quantile(boot$p0, 0.025, na.rm = TRUE),
    p0_high = quantile(boot$p0, 0.975, na.rm = TRUE),
    p1_low = quantile(boot$p1, 0.025, na.rm = TRUE),
    p1_high = quantile(boot$p1, 0.975, na.rm = TRUE),
    rd_low = quantile(boot$rd, 0.025, na.rm = TRUE),
    rd_high = quantile(boot$rd, 0.975, na.rm = TRUE),
    rr_low = quantile(boot$rr, 0.025, na.rm = TRUE),
    rr_high = quantile(boot$rr, 0.975, na.rm = TRUE),
    n_boot = nrow(boot)
  )
}
