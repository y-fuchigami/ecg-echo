source("helpers.R")

root <- "data_raw/mimic-iv-ecg"
machine_path <- find_file(root, "machine_measurements.*\\.csv(\\.gz)?$")
record_path <- find_file(root, "record_list.*\\.csv(\\.gz)?$")

machine <- fread(machine_path, showProgress = FALSE)
machine <- clean_names(machine)

cols <- names(machine)
report_cols <- cols[str_detect(cols, "^report_?[0-9]+$")]
reports <- select(machine, all_of(report_cols))
reports <- as.matrix(reports)

ecg_text <- paste(reports[, 1])
for (report_col in report_cols[-1]) {
  ecg_text <- paste(ecg_text, reports[, report_col])
}
machine$ecg_text <- clean_text(ecg_text)

ecg <- machine
if (!is.na(record_path)) {
  record_list <- fread(record_path, showProgress = FALSE)
  record_list <- clean_names(record_list)

  if ("study_id" %in% names(machine) && "study_id" %in% names(record_list)) {
    ecg <- left_join(machine, record_list, by = "study_id", suffix = c("", "_record"))
  }
}

time_name <- time_col(ecg)
if (!is.na(time_name)) {
  ecg$ecg_time <- to_time(ecg[[time_name]])
} else {
  ecg$ecg_time <- as.POSIXct(NA)
}

ecg <- mutate(
  ecg,
  rad = str_detect(ecg_text, "right axis deviation|rightward axis|right axis"),
  rvh = str_detect(ecg_text, "right ventricular hypertrophy|rvh"),
  rae = str_detect(ecg_text, "right atrial enlargement|right atrial abnormality|right atrial hypertrophy|rae"),
  rv_strain = str_detect(ecg_text, "right ventricular strain|rv strain"),
  pulmonary_pattern = str_detect(ecg_text, "pulmonary disease pattern|pulmonary pattern"),

  rbbb = str_detect(ecg_text, "right bundle branch block|rbbb"),
  incomplete_rbbb = str_detect(ecg_text, "incomplete right bundle branch block|irbbb"),
  s1q3t3 = str_detect(ecg_text, "s1q3t3|s1 q3 t3"),
  anterior_t_abn = str_detect(ecg_text, "anterior t wave abnormality|t wave abnormality.*anterior|anterior.*t wave"),

  right_heart_core = rad | rvh | rae | rv_strain | pulmonary_pattern,
  right_heart_extended = right_heart_core | rbbb | incomplete_rbbb | s1q3t3 | anterior_t_abn
)

fwrite(ecg, "data_processed/ecg.csv", showProgress = FALSE)
