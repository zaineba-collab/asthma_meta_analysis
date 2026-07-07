library(data.table)

# Shared MR configuration helpers.
#
# Supervisor plan:
#   Use already PLINK LD-clumped asthma signals as instruments and do not
#   perform additional clumping in TwoSampleMR.

get_script_context <- function() {
  script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
  project_dir <- dirname(script_dir)

  list(script_file = script_file, script_dir = script_dir, project_dir = project_dir)
}

get_requested_outcome_label <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  label <- if (length(args) >= 1 && nzchar(args[1])) args[1] else Sys.getenv("MR_OUTCOME_LABEL", unset = "")

  if (!nzchar(label)) {
    stop(
      "Outcome label is missing. Use MR_OUTCOME_LABEL=<label> Rscript <script> ",
      "or Rscript <script> <label>.",
      call. = FALSE
    )
  }

  label
}

read_mr_outcome_config <- function(project_dir, outcome_label = get_requested_outcome_label()) {
  requested_label <- outcome_label
  config_file <- file.path(project_dir, "config", "mr_outcomes.tsv")

  if (!file.exists(config_file)) {
    stop("MR outcome config file not found: ", config_file, call. = FALSE)
  }

  config <- fread(config_file)
  required_columns <- c("outcome_label", "outcome_name", "outcome_id", "outcome_type")
  missing_columns <- setdiff(required_columns, names(config))

  if (length(missing_columns) > 0) {
    stop("MR outcome config missing column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  duplicate_labels <- config[duplicated(outcome_label) | duplicated(outcome_label, fromLast = TRUE)]

  if (nrow(duplicate_labels) > 0) {
    stop(
      "Duplicate outcome_label value(s) in MR outcome config: ",
      paste(sort(unique(duplicate_labels$outcome_label)), collapse = ", "),
      call. = FALSE
    )
  }

  selected <- config[config[["outcome_label"]] == requested_label]

  if (nrow(selected) == 0) {
    stop(
      "Outcome label '", requested_label, "' not found in config/mr_outcomes.tsv.",
      call. = FALSE
    )
  }

  if (!nzchar(selected$outcome_id[1]) || identical(selected$outcome_id[1], "TO_FILL")) {
    stop(
      "Outcome ID for label '", requested_label, "' is TO_FILL. Add the real OpenGWAS ID ",
      "to config/mr_outcomes.tsv before running this outcome.",
      call. = FALSE
    )
  }

  as.list(selected[1])
}

print_mr_config <- function(mr_config, output_dir) {
  message("Selected outcome label: ", mr_config$outcome_label)
  message("Outcome name: ", mr_config$outcome_name)
  message("Outcome ID: ", mr_config$outcome_id)
  message("Outcome type: ", mr_config$outcome_type)
  message("Output directory: ", output_dir)
  message("Using already LD-clumped asthma signals. No additional TwoSampleMR clumping performed.")
}
