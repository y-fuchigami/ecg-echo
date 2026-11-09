to01 <- function(x) {
  if (!is.logical(x) && !is.numeric(x)) {
    x <- as.logical(x)
  }
  as.integer(x)
}

primary_outcomes <- c(
  "echo_ph_specific", "echo_ph_primary", "echo_rv_abnormal"
)

trv_outcomes <- c("echo_ph_trv_28", "echo_ph_trv_34")

rv_component_outcomes <- c(
  "rv_dysfunction", "rv_dilation", "ra_enlargement_echo", "tr_moderate_severe"
)

outcomes <- c(primary_outcomes, trv_outcomes, rv_component_outcomes)

outcome_labels <- c(
  echo_ph_specific = "Echo-PH specific",
  echo_ph_primary = "Echo-PH primary",
  echo_rv_abnormal = "Echo RV abnormality",
  echo_ph_trv_28 = "TR velocity >2.8 m/s",
  echo_ph_trv_34 = "TR velocity >3.4 m/s",
  rv_dysfunction = "RV dysfunction",
  rv_dilation = "RV dilation",
  ra_enlargement_echo = "RA enlargement",
  tr_moderate_severe = "Moderate/severe TR"
)

logical_covars <- c(
  "female", "hypertension", "diabetes", "heart_failure", "copd",
  "ckd", "atrial_fibrillation", "cad", "pulmonary_embolism", "icu_at_echo"
)

covars_full_model <- c(
  "age", "female_num", "race_group",
  "hypertension_num", "diabetes_num", "heart_failure_num", "copd_num",
  "ckd_num", "atrial_fibrillation_num", "cad_num",
  "pulmonary_embolism_num", "icu_at_echo_num"
)

read_analysis <- function(path) {
  df <- data.table::fread(path, showProgress = FALSE)
  df$rad_or_rvh_num <- to01(df$rad_or_rvh)
  for (col in c(outcomes, logical_covars)) {
    name <- paste0(col, "_num")
    df[[name]] <- to01(df[[col]])
  }
  df$race_group <- factor(df$race_group)
  df
}
