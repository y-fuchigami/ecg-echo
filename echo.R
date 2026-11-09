source("helpers.R")

ecg_path <- "data_processed/ecg.csv"
measure_path <- "data_raw/mimic-iv-echo/structured-measurement.csv.gz"
study_path <- "data_raw/mimic-iv-echo/echo-study-list.csv"

ecg <- fread(ecg_path, showProgress = FALSE)
ecg <- clean_names(ecg)
echo_raw <- fread(measure_path, showProgress = FALSE)
echo_raw <- clean_names(echo_raw)
echo_study <- fread(study_path, showProgress = FALSE)
echo_study <- clean_names(echo_study)

ecg_keep <- select(
  ecg,
  subject_id,
  study_id,
  ecg_time,
  ecg_text,
  rad,
  rvh,
  rae,
  rv_strain,
  pulmonary_pattern,
  rbbb,
  incomplete_rbbb,
  s1q3t3,
  anterior_t_abn,
  right_heart_core,
  right_heart_extended,
  qrs_axis,
  p_axis,
  t_axis
)
ecg_keep$ecg_time <- to_time(ecg_keep$ecg_time)
ecg_keep <- distinct(ecg_keep)

echo_raw <- select(echo_raw, subject_id, measurement_id, measurement, result, measurement_datetime)
echo_study <- select(echo_study, subject_id, measurement_id, study_id, study_datetime)

echo_raw$measurement <- as.character(echo_raw$measurement)
echo_raw$measurement <- str_squish(echo_raw$measurement)
echo_raw$measurement <- str_to_lower(echo_raw$measurement)

result <- as.character(echo_raw$result)
echo_raw$result_text <- str_squish(result)
echo_raw$result_text <- str_to_lower(echo_raw$result_text)

value <- suppressWarnings(as.numeric(echo_raw$result))
missing <- is.na(value)
text_num <- str_extract(result[missing], "-?\\d+\\.?\\d*")
value[missing] <- suppressWarnings(as.numeric(text_num))
echo_raw$value <- value

echo_raw$measurement_datetime <- to_time(echo_raw$measurement_datetime)
echo_study$study_datetime <- to_time(echo_study$study_datetime)

echo_flags <- mutate(
  echo_raw,
  tr_velocity_value = ifelse(
    measurement %in% "tr_velocity" & !is.na(value) & value >= 0.5 & value <= 6.0,
    value, NA_real_
  ),
  tr_mmhg_value = ifelse(
    measurement %in% "tr_mmhg" & !is.na(value) & value >= 0 & value <= 150,
    value, NA_real_
  ),
  tv_est_pa_press_min_value = ifelse(
    measurement %in% "tv_est_pa_press_min" & !is.na(value) & value >= 0 & value <= 150,
    value, NA_real_
  ),
  lvef_value = ifelse(
    measurement %in% c("lvef", "biplane_lvef", "lvef_3d", "rest_lvef", "rest_biplane_lvef") &
      !is.na(value) & value >= 1 & value <= 100,
    value, NA_real_
  ),
  ivc_diam_value = ifelse(
    measurement %in% "ivc_diam" & !is.na(value) & value >= 0.3 & value <= 5.0,
    value, NA_real_
  ),
  rv_diam_value = ifelse(
    measurement %in% "rv_diam" & !is.na(value) & value >= 0.5 & value <= 8.0,
    value, NA_real_
  ),

  phtn_any_text = measurement %in% c("phtn_severity", "diastolic_phtn", "stress_phtn_severity") &
    str_detect(result_text, "mild|moderate|mod|severe|sev|pulmonary hypertension|phtn"),

  phtn_moderate_severe_text = measurement %in% c("phtn_severity", "diastolic_phtn", "stress_phtn_severity") &
    str_detect(result_text, "moderate|mod|severe|sev"),

  rv_dysfunction = measurement == "rv_function" &
    str_detect(result_text, "mild|moderate|mod|severe|sev|reduced|depressed|hypo|dysfunction|abnormal|poor"),

  rv_modsev_dysfunction = measurement == "rv_function" &
    str_detect(result_text, "moderate|mod|severe|sev|poor"),

  rv_dilation = measurement %in% c("rv_chamber_size", "rv_diam") &
    str_detect(result_text, "mild|moderate|mod|severe|sev|dilat|enlarg|increased"),

  rv_modsev_dilation = measurement %in% c("rv_chamber_size", "rv_diam") &
    str_detect(result_text, "moderate|mod|severe|sev"),

  ra_enlargement_echo = measurement == "ra_size" &
    str_detect(result_text, "mild|moderate|mod|severe|sev|dilat|enlarg|increased"),

  ra_modsev_enlargement_echo = measurement == "ra_size" &
    str_detect(result_text, "moderate|mod|severe|sev"),

  tr_moderate_severe = measurement == "tricuspid_regurg" &
    str_detect(result_text, "moderate|mod|severe|sev"),

  tr_any_abnormal = measurement == "tricuspid_regurg" &
    str_detect(result_text, "trace|mild|moderate|mod|severe|sev"),

  ivc_dilated_text = measurement %in% c("ivc", "ivc_diam") &
    str_detect(result_text, "dilated|dilat|>2.5|> 2.5|>2.1|> 2.1"),

  ivc_dilated_num = !is.na(ivc_diam_value) & ivc_diam_value >= 2.1,

  rv_pressure_overload = str_detect(
    result_text,
    "rv pressure overload|rv pressure/volume overload|right ventricular pressure overload"
  ),

  rv_volume_overload = str_detect(
    result_text,
    "rv pressure/volume overload|rv volume overload|right ventricular volume overload"
  )
)

max_value <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

echo_refined <- group_by(echo_flags, subject_id, measurement_id)
echo_refined <- summarise(
  echo_refined,
  echo_time = if (all(is.na(measurement_datetime))) {
    as.POSIXct(NA, tz = "UTC")
  } else {
    min(measurement_datetime, na.rm = TRUE)
  },

  tr_velocity = max_value(tr_velocity_value),
  tr_mmhg = max_value(tr_mmhg_value),
  tv_est_pa_press_min = max_value(tv_est_pa_press_min_value),
  lvef = max_value(lvef_value),
  ivc_diam = max_value(ivc_diam_value),
  rv_diam = max_value(rv_diam_value),

  phtn_any_text = any(phtn_any_text, na.rm = TRUE),
  phtn_moderate_severe_text = any(phtn_moderate_severe_text, na.rm = TRUE),

  rv_dysfunction = any(rv_dysfunction, na.rm = TRUE),
  rv_modsev_dysfunction = any(rv_modsev_dysfunction, na.rm = TRUE),

  rv_dilation = any(rv_dilation, na.rm = TRUE),
  rv_modsev_dilation = any(rv_modsev_dilation, na.rm = TRUE),

  ra_enlargement_echo = any(ra_enlargement_echo, na.rm = TRUE),
  ra_modsev_enlargement_echo = any(ra_modsev_enlargement_echo, na.rm = TRUE),

  tr_moderate_severe = any(tr_moderate_severe, na.rm = TRUE),
  tr_any_abnormal = any(tr_any_abnormal, na.rm = TRUE),

  ivc_dilated = any(ivc_dilated_text | ivc_dilated_num, na.rm = TRUE),

  rv_pressure_overload = any(rv_pressure_overload, na.rm = TRUE),
  rv_volume_overload = any(rv_volume_overload, na.rm = TRUE),

  n_echo_rows = n(),
  .groups = "drop"
)

echo_study <- distinct(echo_study)

echo_refined <- left_join(echo_refined, echo_study, by = c("subject_id", "measurement_id"))
echo_refined$echo_time_final <- echo_refined$echo_time
missing_time <- is.na(echo_refined$echo_time_final)
echo_refined$echo_time_final[missing_time] <- echo_refined$study_datetime[missing_time]
echo_refined <- mutate(
  echo_refined,
  echo_ph_trv_28 = !is.na(tr_velocity) & tr_velocity > 2.8,
  echo_ph_trv_34 = !is.na(tr_velocity) & tr_velocity > 3.4,
  echo_ph_tr_gradient_31 = !is.na(tr_mmhg) & tr_mmhg >= 31,
  echo_ph_pa_press_40 = !is.na(tv_est_pa_press_min) & tv_est_pa_press_min >= 40,

  echo_ph_primary = echo_ph_trv_28 | echo_ph_tr_gradient_31 | echo_ph_pa_press_40 | phtn_any_text | rv_pressure_overload,

  echo_ph_specific = echo_ph_trv_34 | phtn_moderate_severe_text | rv_pressure_overload,

  echo_rv_abnormal = rv_dysfunction | rv_dilation | ra_enlargement_echo | tr_moderate_severe | ivc_dilated | rv_pressure_overload | rv_volume_overload,

  lvef_reduced = !is.na(lvef) & lvef < 50,
  lvef_preserved = !is.na(lvef) & lvef >= 50
)

ecg_echo_all <- filter(echo_refined, !is.na(echo_time_final))
ecg_keep <- filter(ecg_keep, !is.na(ecg_time))
ecg_echo_all <- inner_join(
  ecg_echo_all,
  ecg_keep,
  by = "subject_id",
  suffix = c("_echo", "_ecg"),
  relationship = "many-to-many"
)
hours <- difftime(ecg_echo_all$echo_time_final, ecg_echo_all$ecg_time, units = "hours")
ecg_echo_all$time_diff_hours <- as.numeric(hours)
ecg_echo_all$abs_time_diff_hours <- abs(ecg_echo_all$time_diff_hours)

ecg_echo_main <- filter(ecg_echo_all, time_diff_hours >= 0, time_diff_hours <= 24 * 7)
ecg_echo_main <- group_by(ecg_echo_main, subject_id, measurement_id)
ecg_echo_main <- slice_min(ecg_echo_main, abs_time_diff_hours, n = 1, with_ties = FALSE)
ecg_echo_main <- ungroup(ecg_echo_main)
fwrite(ecg_echo_main, "data_processed/pre7d.csv", showProgress = FALSE)

for (window_days in c(1, 3, 7)) {
  window_hours <- 24 * window_days
  matched <- filter(ecg_echo_all, abs_time_diff_hours <= window_hours)
  matched <- group_by(matched, subject_id, measurement_id)
  matched <- slice_min(matched, abs_time_diff_hours, n = 1, with_ties = FALSE)
  matched <- ungroup(matched)
  output_path <- paste0("data_processed/pm", window_days, "d.csv")
  fwrite(matched, output_path, showProgress = FALSE)
}
