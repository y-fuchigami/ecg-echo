suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(janitor)
})

find_file <- function(root, pattern) {
  files <- list.files(
    root, pattern = pattern, recursive = TRUE,
    full.names = TRUE, ignore.case = TRUE
  )
  files[1]
}

time_col <- function(df) {
  cols <- c(
    "ecg_time", "study_datetime", "study_time", "charttime",
    "chart_time", "storetime", "echo_time", "datetime",
    "study_date", "date"
  )
  cols <- intersect(cols, names(df))
  cols[1]
}

to_time <- function(x) {
  time <- ymd_hms(x, quiet = TRUE)
  hm <- ymd_hm(x, quiet = TRUE)
  day <- ymd(x, quiet = TRUE)
  day <- as.POSIXct(day, tz = "UTC")

  na <- is.na(time)
  time[na] <- hm[na]
  na <- is.na(time)
  time[na] <- day[na]
  time
}

clean_text <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- str_to_lower(x)
  str_squish(x)
}
