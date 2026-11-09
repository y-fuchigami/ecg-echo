if (!dir.exists("data_processed")) {
  dir.create("data_processed")
}

source("ecg.R", local = new.env())
source("echo.R", local = new.env())
source("clinical.R", local = new.env())
