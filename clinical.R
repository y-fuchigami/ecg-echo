source("helpers.R")

matched_path <- "data_processed/pre7d.csv"

patients_path <- "data_raw/mimic-iv/patients.csv.gz"
admissions_path <- "data_raw/mimic-iv/admissions.csv.gz"
diagnoses_path <- "data_raw/mimic-iv/diagnoses_icd.csv.gz"
icustays_path <- "data_raw/mimic-iv/icustays.csv.gz"

df <- fread(matched_path, showProgress = FALSE)
df <- clean_names(df)

patients <- fread(patients_path, showProgress = FALSE)
patients <- clean_names(patients)

admissions <- fread(admissions_path, showProgress = FALSE)
admissions <- clean_names(admissions)

diagnoses <- fread(diagnoses_path, showProgress = FALSE)
diagnoses <- clean_names(diagnoses)

icustays <- fread(icustays_path, showProgress = FALSE)
icustays <- clean_names(icustays)

df <- mutate(
  df,
  echo_time_final = to_time(echo_time_final),
  ecg_time = to_time(ecg_time)
)

patients <- mutate(patients, dod = to_time(dod))

admissions <- mutate(
  admissions,
  admittime = to_time(admittime),
  dischtime = to_time(dischtime),
  deathtime = to_time(deathtime)
)

icustays <- mutate(
  icustays,
  intime = to_time(intime),
  outtime = to_time(outtime)
)

df <- mutate(
  df,
  qrs_axis_num = suppressWarnings(as.numeric(qrs_axis)),
  qrs_axis_right_90 = !is.na(qrs_axis_num) & qrs_axis_num > 90,
  qrs_axis_right_110 = !is.na(qrs_axis_num) & qrs_axis_num > 110,

  right_heart_core_no_rae = rad | rvh | rv_strain | pulmonary_pattern,
  rad_or_rvh = rad | rvh,
  isolated_rbbb = rbbb & !right_heart_core,

  right_heart_pattern = right_heart_core_no_rae
)

admitted <- left_join(df, admissions, by = "subject_id", relationship = "many-to-many")
admitted <- filter(admitted, echo_time_final >= admittime, echo_time_final <= dischtime)
hours <- difftime(admitted$echo_time_final, admitted$admittime, units = "hours")
admitted$admit_hours <- as.numeric(hours)
admitted <- group_by(admitted, subject_id, measurement_id)
admitted <- slice_min(admitted, admit_hours, n = 1, with_ties = FALSE)
admitted <- ungroup(admitted)

adm_ids <- select(admitted, subject_id, measurement_id)
no_adm <- anti_join(df, adm_ids, by = c("subject_id", "measurement_id"))
no_adm <- mutate(
  no_adm,
  hadm_id = NA_integer_,
  admittime = as.POSIXct(NA),
  dischtime = as.POSIXct(NA),
  deathtime = as.POSIXct(NA),
  admission_type = NA_character_,
  admission_location = NA_character_,
  discharge_location = NA_character_,
  insurance = NA_character_,
  language = NA_character_,
  marital_status = NA_character_,
  race = NA_character_,
  hospital_expire_flag = NA_integer_
)

cols <- intersect(names(admitted), names(no_adm))
admitted <- select(admitted, all_of(cols))
no_adm <- select(no_adm, all_of(cols))
df <- bind_rows(admitted, no_adm)

patients <- select(patients, subject_id, gender, anchor_age, anchor_year, dod)
df <- left_join(df, patients, by = "subject_id")
df <- mutate(
  df,
  echo_year = year(echo_time_final),
  age = anchor_age + (echo_year - anchor_year),
  age = ifelse(age < 18 | age > 120, NA_real_, age),
  female = gender == "F"
)

race_group <- character(nrow(df))
for (i in seq_len(nrow(df))) {
  race <- str_to_lower(df$race[i])
  if (is.na(race)) {
    group <- "Unknown"
  } else if (str_detect(race, "white")) {
    group <- "White"
  } else if (str_detect(race, "black|african")) {
    group <- "Black"
  } else if (str_detect(race, "asian")) {
    group <- "Asian"
  } else if (str_detect(race, "hispanic|latino")) {
    group <- "Hispanic"
  } else {
    group <- "Other"
  }
  race_group[i] <- group
}
df$race_group <- race_group

icu <- select(df, subject_id, measurement_id, echo_time_final)
icu <- left_join(icu, icustays, by = "subject_id", relationship = "many-to-many")
icu <- mutate(icu, icu_at_echo = echo_time_final >= intime & echo_time_final <= outtime)
icu <- group_by(icu, subject_id, measurement_id)
icu <- summarise(
  icu,
  icu_at_echo = any(icu_at_echo, na.rm = TRUE),
  .groups = "drop"
)

df <- left_join(df, icu, by = c("subject_id", "measurement_id"))

diagnoses <- mutate(
  diagnoses,
  code = str_to_upper(icd_code),
  code = str_replace_all(code, "\\.", "")
)

dx <- mutate(
  diagnoses,
  hypertension = ifelse(
    icd_version == 9,
    str_detect(code, "^401|^402|^403|^404|^405"),
    str_detect(code, "^I10|^I11|^I12|^I13|^I15")
  ),

  diabetes = ifelse(
    icd_version == 9,
    str_detect(code, "^250"),
    str_detect(code, "^E10|^E11|^E13|^E14")
  ),

  heart_failure = ifelse(
    icd_version == 9,
    str_detect(code, "^428"),
    str_detect(code, "^I50")
  ),

  copd = ifelse(
    icd_version == 9,
    str_detect(code, "^490|^491|^492|^496"),
    str_detect(code, "^J40|^J41|^J42|^J43|^J44")
  ),

  ckd = ifelse(
    icd_version == 9,
    str_detect(code, "^585"),
    str_detect(code, "^N18")
  ),

  atrial_fibrillation = ifelse(
    icd_version == 9,
    str_detect(code, "^42731"),
    str_detect(code, "^I48")
  ),

  cad = ifelse(
    icd_version == 9,
    str_detect(code, "^410|^411|^412|^413|^414"),
    str_detect(code, "^I20|^I21|^I22|^I23|^I24|^I25")
  ),

  pulmonary_embolism = ifelse(
    icd_version == 9,
    str_detect(code, "^4151"),
    str_detect(code, "^I26")
  )
)
dx <- group_by(dx, subject_id, hadm_id)
dx <- summarise(
  dx,
  hypertension = any(hypertension, na.rm = TRUE),
  diabetes = any(diabetes, na.rm = TRUE),
  heart_failure = any(heart_failure, na.rm = TRUE),
  copd = any(copd, na.rm = TRUE),
  ckd = any(ckd, na.rm = TRUE),
  atrial_fibrillation = any(atrial_fibrillation, na.rm = TRUE),
  cad = any(cad, na.rm = TRUE),
  pulmonary_embolism = any(pulmonary_embolism, na.rm = TRUE),
  .groups = "drop"
)

df <- left_join(df, dx, by = c("subject_id", "hadm_id"))

dx_names <- c(
  "hypertension", "diabetes", "heart_failure", "copd",
  "ckd", "atrial_fibrillation", "cad", "pulmonary_embolism"
)
for (name in dx_names) {
  values <- df[[name]]
  values[is.na(values)] <- FALSE
  df[[name]] <- values
}

fwrite(df, "data_processed/analysis.csv", showProgress = FALSE)
