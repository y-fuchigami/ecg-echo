library(data.table)
source("data.R")
source("clustered.R")
source("diagnostics.R")

if (!dir.exists(output)) {
  dir.create(output)
}
analysis_path <- file.path(
  project, "data_processed/analysis.csv"
)
df <- read_analysis(analysis_path)

save_table <- function(table, file) {
  fwrite(table, file.path(output, file), showProgress = FALSE)
}

number <- function(value) {
  format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
}

count_pct <- function(count, total) {
  sprintf("%s (%.1f%%)", number(count), 100 * count / total)
}

percent_ci <- function(estimate, lower, upper) {
  sprintf("%.1f (%.1f–%.1f)", 100 * estimate, 100 * lower, 100 * upper)
}

ratio_ci <- function(estimate, lower, upper) {
  sprintf("%.2f (%.2f–%.2f)", estimate, lower, upper)
}

variables <- c(
  "N" = "n", "Age, median [IQR]" = "age", "Female" = "female",
  "ICU at echo" = "icu_at_echo", "Hypertension" = "hypertension",
  "Diabetes" = "diabetes", "Heart failure" = "heart_failure", "COPD" = "copd",
  "Chronic kidney disease" = "ckd", "Atrial fibrillation" = "atrial_fibrillation",
  "Coronary artery disease" = "cad", "Pulmonary embolism" = "pulmonary_embolism",
  "Echo-PH primary" = "echo_ph_primary", "Echo-PH specific" = "echo_ph_specific",
  "Echo RV abnormality" = "echo_rv_abnormal", "RV dysfunction" = "rv_dysfunction",
  "RV dilation" = "rv_dilation", "RA enlargement" = "ra_enlargement_echo",
  "Moderate/severe TR" = "tr_moderate_severe"
)
baseline <- function(data) {
  result <- character(length(variables))
  names(result) <- names(variables)
  for (i in seq_along(variables)) {
    variable <- variables[i]
    if (variable == "n") {
      result[i] <- number(nrow(data))
    } else if (variable == "age") {
      med <- median(data$age, na.rm = TRUE)
      lo <- quantile(data$age, 0.25, na.rm = TRUE)
      hi <- quantile(data$age, 0.75, na.rm = TRUE)
      result[i] <- sprintf(
        "%.0f [%.0f–%.0f]", med, lo, hi
      )
    } else {
      count <- sum(data[[variable]], na.rm = TRUE)
      result[i] <- count_pct(count, nrow(data))
    }
  }
  result
}
keep <- df$rad_or_rvh_num == 0
negative <- df[keep]
keep <- df$rad_or_rvh_num == 1
positive <- df[keep]
table1 <- data.table(
  Characteristic = names(variables), Overall = baseline(df),
  negative = baseline(negative), positive = baseline(positive)
)
save_table(table1, "Table1.csv")

exposures <- c(
  rad_or_rvh = "RAD/RVH", rad = "RAD", rvh = "RVH",
  qrs_axis_right_110 = "QRS axis >110°", rbbb = "RBBB"
)
diagnostic_cells <- function(diagnostics, component = FALSE) {
  proportion_metrics <- c("sensitivity", "specificity", "ppv", "npv")
  ratio_metrics <- c("lr_pos", "lr_neg")
  metrics <- c(proportion_metrics, ratio_metrics)
  result <- data.table(
    n = number(diagnostics$n), tp = number(diagnostics$tp),
    fp = number(diagnostics$fp), fn = number(diagnostics$fn),
    tn = number(diagnostics$tn)
  )
  for (metric in metrics) {
    estimate <- diagnostics[[metric]]
    lower <- diagnostics[[paste0(metric, "_low")]]
    upper <- diagnostics[[paste0(metric, "_high")]]
    if (metric %in% proportion_metrics) {
      result[[metric]] <- percent_ci(estimate, lower, upper)
    } else {
      result[[metric]] <- ratio_ci(estimate, lower, upper)
    }
  }
  if (component) {
    result <- cbind(Outcome = diagnostics$outcome, result)
  } else {
    result <- cbind(
      ECG = diagnostics$exposure, Outcome = diagnostics$outcome, result
    )
  }
  result
}
diagnostics <- list()
for (outcome in primary_outcomes) {
  for (exposure in names(exposures)) {
    diagnostics[[length(diagnostics) + 1L]] <- diagnostic_one(
      df, exposure, outcome, exposures[[exposure]],
      outcome_labels[[outcome]], ci_method
    )
  }
}
diagnostics <- rbindlist(diagnostics)
table2 <- diagnostic_cells(diagnostics)
save_table(table2, "Table2.csv")

clustered_models <- list()
for (outcome in outcomes) {
  clustered_models[[outcome]] <- fit_clustered_logit(df, outcome)
}
clustered_models <- rbindlist(clustered_models)
p_text <- function(p) {
  ifelse(p < .001, "<0.001", sprintf("%.3f", p))
}
table3 <- data.table(
  Outcome = outcome_labels[clustered_models$outcome],
  Pairs = number(clustered_models$n_pairs),
  Patients = number(clustered_models$n_patients),
  Events = number(clustered_models$n_events),
  Positive = number(clustered_models$n_exposed),
  OR = ratio_ci(
    clustered_models$OR, clustered_models$CI_low, clustered_models$CI_high
  ),
  P = p_text(clustered_models$p_value)
)
save_table(table3, "Table3.csv")

df$echo_time_final <- as.POSIXct(df$echo_time_final, tz = "UTC")
df$ecg_time <- as.POSIXct(df$ecg_time, tz = "UTC")
time_interval <- difftime(
  df$echo_time_final, df$ecg_time, units = "hours"
)
df$absolute_interval <- abs(as.numeric(time_interval))
right_heart_num <- function(df) {
  pattern <- df$rad | df$rvh | df$rv_strain | df$pulmonary_pattern
  as.integer(pattern)
}
df$right_heart_pattern_num <- right_heart_num(df)
for (col in c(logical_covars, primary_outcomes, trv_outcomes)) {
  name <- paste0(col, "_01")
  df[[name]] <- as.integer(df[[col]] == TRUE)
}

model_labels <- c(
  "Model 1: unadjusted", "Model 2: age + sex + race",
  "Model 3: demographics + comorbidities + ICU"
)
model_values <- function(fit) {
  interval <- suppressMessages(confint(fit, parm = "right_heart_pattern_num"))
  coefs <- coef(summary(fit))
  estimate <- coefs["right_heart_pattern_num", "Estimate"]
  p <- coefs["right_heart_pattern_num", "Pr(>|z|)"]
  lower <- exp(interval[["2.5 %"]])
  upper <- exp(interval[["97.5 %"]])
  list(OR = ratio_ci(exp(estimate), lower, upper), P = p_text(p))
}

fit_three <- function(data, col) {
  model_terms <- c(
    "right_heart_pattern_num",
    "right_heart_pattern_num + age + female_01 + race_group",
    paste(
      "right_heart_pattern_num + age + female_01 + race_group +",
      "hypertension_01 + diabetes_01 + heart_failure_01 + copd_01 + ckd_01 +",
      "atrial_fibrillation_01 + cad_01 + pulmonary_embolism_01 + icu_at_echo_01"
    )
  )
  models <- list()
  for (i in seq_along(model_terms)) {
    form <- paste(col, "~", model_terms[i])
    form <- as.formula(form)
    fit <- glm(form, data = data, family = binomial())
    values <- model_values(fit)
    models[[i]] <- data.table(
      Model = model_labels[i], OR = values$OR,
      P = values$P, Model_pairs = number(nobs(fit))
    )
  }
  rbindlist(models)
}

rows <- order(df$subject_id, df$echo_time_final)
first <- df[rows]
keep <- !duplicated(first$subject_id)
first <- first[keep]
rows <- order(df$subject_id, df$absolute_interval)
closest <- df[rows]
keep <- !duplicated(closest$subject_id)
closest <- closest[keep]
patient_sets <- list(
  "Closest ECG–echo per patient" = closest, "First echo per patient" = first
)
s1 <- list()
for (label in names(patient_sets)) {
  patient_data <- patient_sets[[label]]
  for (outcome in primary_outcomes) {
    col <- paste0(outcome, "_01")
    models <- fit_three(patient_data, col)
    s1[[length(s1) + 1L]] <- data.table(
      Dataset = label, Outcome = outcome_labels[[outcome]], models,
      Input_patients = number(uniqueN(patient_data$subject_id)),
      Input_positive = number(sum(patient_data$right_heart_pattern_num == 1)),
      Input_events = number(sum(patient_data[[col]] == 1))
    )
  }
}
s1 <- rbindlist(s1)
save_table(s1, "S1.csv")

windows <- c(
  "Pre-Echo ECG within 7 days" = "pre7d", "ECG within +/- 1 day" = "pm1d",
  "ECG within +/- 3 days" = "pm3d", "ECG within +/- 7 days" = "pm7d"
)
s2 <- list()
for (label in names(windows)) {
  window_file <- paste0(windows[[label]], ".csv")
  window_path <- file.path(project, "data_processed", window_file)
  window_data <- fread(window_path, showProgress = FALSE)
  window_data$right_heart_pattern_num <- right_heart_num(window_data)
  for (outcome in primary_outcomes) {
    window_data$outcome_num <- as.integer(window_data[[outcome]] == TRUE)
    metric <- diagnostic_one(window_data, "right_heart_pattern_num", "outcome_num")
    fit <- glm(
      outcome_num ~ right_heart_pattern_num, data = window_data, family = binomial()
    )
    values <- model_values(fit)
    s2[[length(s2) + 1L]] <- data.table(
      Dataset = label, Outcome = outcome_labels[[outcome]],
      Input_pairs = number(nrow(window_data)),
      Input_patients = number(uniqueN(window_data$subject_id)),
      Input_positive = number(sum(window_data$right_heart_pattern_num == 1)),
      Input_events = number(sum(window_data$outcome_num == 1)),
      Sensitivity = sprintf("%.1f", 100 * metric$sensitivity),
      Specificity = sprintf("%.1f", 100 * metric$specificity),
      PPV = sprintf("%.1f", 100 * metric$ppv), NPV = sprintf("%.1f", 100 * metric$npv),
      OR = values$OR, P = values$P
    )
  }
}
s2 <- rbindlist(s2)
save_table(s2, "S2.csv")

keep <- !is.na(df$tr_velocity)
trv_data <- df[keep]
s3a <- list()
s3b <- list()
for (outcome in trv_outcomes) {
  col <- paste0(outcome, "_01")
  metric <- diagnostic_one(trv_data, "right_heart_pattern_num", col)
  s3a[[outcome]] <- data.table(
    Outcome = outcome_labels[[outcome]], N = number(nrow(trv_data)),
    Patients = number(uniqueN(trv_data$subject_id)),
    Positive = number(sum(trv_data$right_heart_pattern_num == 1)),
    Events = number(sum(trv_data[[col]] == 1)),
    Sensitivity = sprintf("%.1f", 100 * metric$sensitivity),
    Specificity = sprintf("%.1f", 100 * metric$specificity),
    PPV = sprintf("%.1f", 100 * metric$ppv), NPV = sprintf("%.1f", 100 * metric$npv),
    LRpos = sprintf("%.2f", metric$lr_pos), LRneg = sprintf("%.2f", metric$lr_neg)
  )
  s3b[[outcome]] <- cbind(
    Outcome = outcome_labels[[outcome]], fit_three(trv_data, col)
  )
}
s3a <- rbindlist(s3a)
save_table(s3a, "S3a.csv")
s3b <- rbindlist(s3b)
save_table(s3b, "S3b.csv")

marginal <- fread(bootstrap_path, showProgress = FALSE)
row_order <- match(outcomes, marginal$outcome)
marginal <- marginal[row_order]
s4 <- data.table(
  Outcome = outcome_labels[marginal$outcome],
  p0 = percent_ci(marginal$adjusted_probability_unexposed, marginal$p0_low, marginal$p0_high),
  p1 = percent_ci(marginal$adjusted_probability_exposed, marginal$p1_low, marginal$p1_high),
  rd = percent_ci(marginal$adjusted_risk_difference, marginal$rd_low, marginal$rd_high),
  rr = ratio_ci(marginal$adjusted_risk_ratio, marginal$rr_low, marginal$rr_high),
  Pairs = number(marginal$n_pairs), Patients = number(marginal$n_patients),
  Events = number(marginal$n_events), Bootstrap = marginal$n_boot
)
save_table(s4, "S4.csv")

s5 <- list()
for (outcome in rv_component_outcomes) {
  s5[[outcome]] <- diagnostic_one(
    df, "rad_or_rvh", outcome,
    label_outcome = outcome_labels[[outcome]], ci_method = ci_method
  )
}
s5 <- rbindlist(s5)
s5 <- diagnostic_cells(s5, component = TRUE)
save_table(s5, "S5.csv")

review <- fread(review_path, select = c(
  "validation_id", "automated_rad", "automated_rvh", "automated_rad_or_rvh",
  "manual_rad_present", "manual_rvh_present", "manual_rad_or_rvh_present",
  "manual_uncertain_or_borderline_only", "manual_negated_only"
), showProgress = FALSE)
targets <- c("RAD/RVH composite" = "rad_or_rvh", "RAD" = "rad", "RVH" = "rvh")
s6 <- list()
for (label in names(targets)) {
  target <- targets[[label]]
  automated <- to01(review[[paste0("automated_", target)]])
  manual <- to01(review[[paste0("manual_", target, "_present")]])
  tp <- sum(automated == 1 & manual == 1)
  fp <- sum(automated == 1 & manual == 0)
  fn <- sum(automated == 0 & manual == 1)
  tn <- sum(automated == 0 & manual == 0)
  numerators <- c(
    PPV = tp, NPV = tn, Sensitivity = tp,
    Specificity = tn, Accuracy = tp + tn
  )
  denominators <- c(
    PPV = tp + fp, NPV = tn + fn,
    Sensitivity = tp + fn,
    Specificity = tn + fp, Accuracy = length(automated)
  )
  metrics <- character(length(numerators))
  names(metrics) <- names(numerators)
  for (metric in names(numerators)) {
    numerator <- numerators[[metric]]
    denominator <- denominators[[metric]]
    interval <- binom.test(numerator, denominator)$conf.int
    metrics[[metric]] <- percent_ci(
      numerator / denominator, interval[1], interval[2]
    )
  }
  s6[[label]] <- data.table(
    Target = label, Reviewed = number(length(automated)),
    TP = number(tp), FP = number(fp),
    FN = number(fn), TN = number(tn),
    PPV = metrics[["PPV"]], NPV = metrics[["NPV"]],
    Sensitivity = metrics[["Sensitivity"]], Specificity = metrics[["Specificity"]],
    Accuracy = metrics[["Accuracy"]]
  )
}
s6 <- rbindlist(s6)
save_table(s6, "S6.csv")
different <- review$automated_rad_or_rvh != review$manual_rad_or_rvh_present
discordant <- review[different]
discordant <- discordant[order(discordant$validation_id, na.last = FALSE)]
save_table(discordant, "S7.csv")

absolute_interval <- df$absolute_interval
interval_median <- median(absolute_interval)
interval_lower <- quantile(absolute_interval, .25)
interval_upper <- quantile(absolute_interval, .75)
s8 <- data.table(
  Item = c(
    "ECG–echocardiography pairs", "Patients", "RAD/RVH-positive pairs",
    "TR velocity available", "TR velocity missing",
    "Median absolute ECG–echocardiography interval, hours",
    "Interquartile range of absolute interval, hours",
    "ECG–echocardiography pairs within 24 hours",
    "ECG–echocardiography pairs within 6 hours"
  ),
  Value = c(
    number(nrow(df)), number(uniqueN(df$subject_id)),
    number(sum(df$rad_or_rvh_num)),
    count_pct(sum(!is.na(df$tr_velocity)), nrow(df)),
    count_pct(sum(is.na(df$tr_velocity)), nrow(df)),
    sprintf("%.1f", interval_median),
    sprintf("%.1f–%.1f", interval_lower, interval_upper),
    count_pct(sum(absolute_interval <= 24), nrow(df)),
    count_pct(sum(absolute_interval <= 6), nrow(df))
  )
)
save_table(s8, "S8.csv")
