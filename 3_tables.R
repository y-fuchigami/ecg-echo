project <- "."
output <- "tables"
bootstrap_path <- "results/marginal.csv"
review_path <- file.path(project, "output/revision1/manual_validation/manual_validation_completed_merged.csv")
ci_method <- "exact"

source("tables.R", local = TRUE)
