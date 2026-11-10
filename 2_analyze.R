source("data.R")
source("clustered.R")
source("diagnostics.R")
source("marginal.R")

input <- "data_processed/analysis.csv"
output <- "results"
B <- 300
ci_method <- "wilson"

if (!dir.exists(output)) {
  dir.create(output)
}
df <- read_analysis(input)

clustered <- list()
for (outcome in outcomes) {
  clustered[[outcome]] <- fit_clustered_logit(df, outcome)
}
clustered <- data.table::rbindlist(clustered)
data.table::fwrite(clustered, file.path(output, "clustered.csv"), showProgress = FALSE)

exposures <- c(
  rad_or_rvh = "Machine-reported RAD/RVH",
  rad = "Machine-reported RAD",
  rvh = "Machine-reported RVH",
  qrs_axis_right_110 = "QRS axis >110°",
  rbbb = "Machine-reported RBBB"
)

diagnostics <- list()
for (exposure in names(exposures)) {
  for (outcome in primary_outcomes) {
    diagnostics[[length(diagnostics) + 1L]] <- diagnostic_one(
      df, exposure, outcome, exposures[[exposure]],
      outcome_labels[[outcome]], ci_method
    )
  }
}
diagnostics <- data.table::rbindlist(diagnostics)
data.table::fwrite(diagnostics, file.path(output, "diagnostics.csv"), showProgress = FALSE)

marginal <- list()
for (outcome in outcomes) {
  marginal[[outcome]] <- fit_marginal_probs(df, outcome, B)
}
marginal <- data.table::rbindlist(marginal)
data.table::fwrite(marginal, file.path(output, "marginal.csv"), showProgress = FALSE)
