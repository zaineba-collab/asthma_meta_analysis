#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# Mode estimators use bootstrap sampling; fix the RNG seed for reproducibility.
set.seed(2401)

# Primary MR only: Asthma -> 12 local FinnGen outcomes.
# Inputs are the saved MR-ready datasets; this script does not harmonise,
# clump, query external services, run sensitivity analyses, or create plots.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(arg)) sub("^--file=", "", arg[1]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)

ready_dir <- p("results", "mr", "finngen", "harmonised", "mr_ready")
config_file <- p("config", "finngen_mr_outcomes.tsv")
harm_qc_file <- p("results", "mr", "finngen", "qc", "finngen_harmonisation_qc.tsv")
result_dir <- p("results", "mr", "finngen", "results", "primary")
qc_dir <- p("results", "mr", "finngen", "qc")
error_file <- file.path(qc_dir, "finngen_primary_mr_errors.tsv")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(config_file)) stop("Missing outcome config: ", config_file, call. = FALSE)
if (!file.exists(harm_qc_file)) stop("Missing harmonisation QC: ", harm_qc_file, call. = FALSE)

cfg <- fread(config_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
harm_qc <- fread(harm_qc_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
if (nrow(cfg) != 12L || uniqueN(cfg$outcome_id) != 12L) stop("Outcome config must contain 12 unique outcomes", call. = FALSE)

methods <- data.table(
  method_id = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"),
  method = c("Inverse variance weighted", "MR Egger", "Weighted median", "Simple mode", "Weighted mode"),
  method_order = 1:5
)
installed_methods <- as.data.table(TwoSampleMR::mr_method_list())
missing_methods <- setdiff(methods$method_id, installed_methods$obj)
if (length(missing_methods)) stop("Requested TwoSampleMR methods unavailable: ", paste(missing_methods, collapse = ", "), call. = FALSE)

required <- c("SNP", "mr_keep", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome",
              "effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome",
              "other_allele.outcome", "exposure", "outcome")
errors <- list()
results <- list()
add_error <- function(outcome, outcome_id, stage, method_id = NA_character_, message) {
  errors[[length(errors) + 1L]] <<- data.table(
    outcome = as.character(outcome), outcome_id = as.character(outcome_id),
    stage = stage, method_id = method_id, error = as.character(message)
  )
}
failed_rows <- function(meta, message) {
  methods[, .(
    exposure = "Asthma", outcome = as.character(meta$outcome), outcome_id = as.character(meta$outcome_id),
    outcome_type = as.character(meta$outcome_type), method, method_id, nsnp = NA_integer_,
    beta = NA_real_, SE = NA_real_, p_value = NA_real_, beta_lower_95CI = NA_real_,
    beta_upper_95CI = NA_real_, OR = NA_real_, OR_lower_95CI = NA_real_, OR_upper_95CI = NA_real_,
    qc_pass = FALSE, qc_message = message, method_order
  )]
}

for (i in seq_len(nrow(cfg))) {
  meta <- cfg[i]
  input_file <- file.path(ready_dir, paste0(meta$output_stub, "_mr_ready.tsv"))
  message("Primary MR: Asthma -> ", meta$outcome)

  outcome_result <- tryCatch({
    if (!file.exists(input_file)) stop("Missing MR-ready input: ", input_file)
    dat <- fread(input_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
    absent <- setdiff(required, names(dat))
    if (length(absent)) stop("Missing fields: ", paste(absent, collapse = ", "))
    if (anyDuplicated(dat$SNP)) stop("SNP is not unique")
    if (nrow(dat) == 0L || any(is.na(dat$mr_keep)) || any(!dat$mr_keep)) stop("Not all rows have mr_keep == TRUE")
    numeric_fields <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
    if (!all(vapply(dat[, ..numeric_fields], is.numeric, logical(1)))) stop("Beta/SE fields are not numeric")
    if (anyNA(dat[, ..numeric_fields])) stop("Missing beta or SE")
    if (any(!is.finite(dat$beta.exposure)) || any(!is.finite(dat$beta.outcome))) stop("Non-finite beta")
    if (any(!is.finite(dat$se.exposure)) || any(!is.finite(dat$se.outcome)) ||
        any(dat$se.exposure <= 0) || any(dat$se.outcome <= 0)) stop("SE is non-finite or non-positive")
    allele_fields <- c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome")
    if (anyNA(dat[, ..allele_fields]) || any(dat$effect_allele.exposure != dat$effect_allele.outcome) ||
        any(dat$other_allele.exposure != dat$other_allele.outcome)) stop("Exposure and outcome alleles are not harmonised")
    expected <- harm_qc[outcome_id == as.character(meta$outcome_id), final_mr_keep_snps]
    if (length(expected) != 1L) stop("No unique harmonisation QC count for outcome")
    if (nrow(dat) != expected) stop("MR-ready row count ", nrow(dat), " differs from harmonisation QC count ", expected)

    method_rows <- vector("list", nrow(methods))
    for (j in seq_len(nrow(methods))) {
      method_meta <- methods[j]
      method_rows[[j]] <- tryCatch({
        fit <- as.data.table(TwoSampleMR::mr(as.data.frame(dat), method_list = method_meta$method_id))
        if (nrow(fit) != 1L) stop("Estimator returned ", nrow(fit), " rows")
        b <- as.numeric(fit$b[1]); se <- as.numeric(fit$se[1]); pv <- as.numeric(fit$pval[1]); n <- as.integer(fit$nsnp[1])
        lo <- b - 1.96 * se; hi <- b + 1.96 * se
        is_binary <- identical(as.character(meta$outcome_type), "binary")
        or <- if (is_binary) exp(b) else NA_real_
        or_lo <- if (is_binary) exp(lo) else NA_real_
        or_hi <- if (is_binary) exp(hi) else NA_real_
        checks <- c(
          !is.na(n) && n > 0L && n <= nrow(dat),
          is.finite(b), is.finite(se) && se > 0,
          is.finite(pv) && pv >= 0 && pv <= 1,
          is.finite(lo) && is.finite(hi) && lo <= b && b <= hi,
          !is_binary || (is.finite(or) && or > 0 && is.finite(or_lo) && is.finite(or_hi) && or_lo <= or && or <= or_hi)
        )
        labels <- c("invalid nsnp", "non-finite beta", "invalid SE", "invalid p-value", "invalid beta CI", "invalid OR/OR CI")
        qc_message <- if (all(checks)) "PASS" else paste(labels[!checks], collapse = "; ")
        if (!all(checks)) add_error(meta$outcome, meta$outcome_id, "result_qc", method_meta$method_id, qc_message)
        if (!is.na(n) && n != nrow(dat)) {
          nsnp_message <- paste0("nsnp ", n, " differs from MR-ready input count ", nrow(dat))
          add_error(meta$outcome, meta$outcome_id, "nsnp_qc", method_meta$method_id, nsnp_message)
          qc_message <- if (identical(qc_message, "PASS")) nsnp_message else paste(qc_message, nsnp_message, sep = "; ")
        }
        data.table(
          exposure = "Asthma", outcome = as.character(meta$outcome), outcome_id = as.character(meta$outcome_id),
          outcome_type = as.character(meta$outcome_type), method = method_meta$method,
          method_id = method_meta$method_id, nsnp = n, beta = b, SE = se, p_value = pv,
          beta_lower_95CI = lo, beta_upper_95CI = hi, OR = or,
          OR_lower_95CI = or_lo, OR_upper_95CI = or_hi,
          qc_pass = all(checks) && n == nrow(dat), qc_message = qc_message,
          method_order = method_meta$method_order
        )
      }, error = function(e) {
        msg <- conditionMessage(e)
        add_error(meta$outcome, meta$outcome_id, "estimator", method_meta$method_id, msg)
        failed_rows(meta, msg)[method_id == method_meta$method_id]
      })
    }
    rbindlist(method_rows, fill = TRUE)
  }, error = function(e) {
    msg <- conditionMessage(e)
    add_error(meta$outcome, meta$outcome_id, "input_validation", message = msg)
    failed_rows(meta, msg)
  })

  setorder(outcome_result, method_order)
  fwrite(outcome_result[, !"method_order"], file.path(result_dir, paste0(meta$output_stub, "_primary_mr.tsv")),
         sep = "\t", quote = FALSE, na = "NA")
  results[[i]] <- outcome_result
}

all_results <- rbindlist(results, fill = TRUE)
all_results[, outcome_order := match(outcome_id, cfg$outcome_id)]
setorder(all_results, outcome_order, method_order)
fwrite(all_results[, !c("outcome_order", "method_order")],
       file.path(result_dir, "finngen_primary_mr_all_methods.tsv"), sep = "\t", quote = FALSE, na = "NA")

ivw <- copy(all_results[method_id == "mr_ivw"])
bonferroni_threshold <- 0.05 / 12
ivw[, `:=`(
  p_bonferroni = pmin(p_value * 12, 1),
  bonferroni_significant = !is.na(p_value) & p_value < bonferroni_threshold,
  p_fdr_bh = p.adjust(p_value, method = "BH")
)]
ivw[, fdr_significant := !is.na(p_fdr_bh) & p_fdr_bh < 0.05]
ivw_summary <- ivw[, .(outcome, outcome_id, outcome_type, nsnp, beta, SE, beta_lower_95CI,
  beta_upper_95CI, p_value, OR, OR_lower_95CI, OR_upper_95CI, p_bonferroni,
  bonferroni_significant, p_fdr_bh, fdr_significant)]
fwrite(ivw_summary, file.path(result_dir, "finngen_primary_ivw_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

direction <- function(x) fifelse(is.na(x), "NA", fifelse(x > 0, "positive", fifelse(x < 0, "negative", "zero")))
wide <- dcast(all_results, outcome + outcome_id ~ method_id, value.var = "beta")
wide_p <- dcast(all_results, outcome + outcome_id ~ method_id, value.var = "p_value")
consistency <- wide[, .(
  outcome, outcome_id,
  ivw_beta_direction = direction(mr_ivw),
  mr_egger_beta_direction = direction(mr_egger_regression),
  weighted_median_beta_direction = direction(mr_weighted_median),
  simple_mode_beta_direction = direction(mr_simple_mode),
  weighted_mode_beta_direction = direction(mr_weighted_mode),
  methods_positive = rowSums(.SD > 0, na.rm = TRUE),
  methods_negative = rowSums(.SD < 0, na.rm = TRUE),
  all_five_same_direction = rowSums(!is.na(.SD)) == 5L & (rowSums(.SD > 0, na.rm = TRUE) == 5L | rowSums(.SD < 0, na.rm = TRUE) == 5L)
), .SDcols = methods$method_id]
consistency <- merge(consistency, wide_p[, .(outcome_id, ivw_p_value = mr_ivw,
  weighted_median_p_value = mr_weighted_median, mr_egger_p_value = mr_egger_regression)], by = "outcome_id", sort = FALSE)
consistency[, outcome_order := match(outcome_id, cfg$outcome_id)]
setorder(consistency, outcome_order)
setcolorder(consistency, c("outcome", "outcome_id", setdiff(names(consistency), c("outcome", "outcome_id", "outcome_order")), "outcome_order"))
fwrite(consistency[, !"outcome_order"], file.path(result_dir, "finngen_mr_method_consistency.tsv"), sep = "\t", quote = FALSE, na = "NA")

if (length(errors)) {
  fwrite(rbindlist(errors, fill = TRUE), error_file, sep = "\t", quote = FALSE, na = "NA")
} else {
  fwrite(data.table(outcome = character(), outcome_id = character(), stage = character(),
                    method_id = character(), error = character()),
         error_file, sep = "\t", quote = FALSE, na = "NA")
}

message("Bonferroni threshold (0.05/12): ", format(bonferroni_threshold, digits = 16))
message("Primary method estimates passing QC: ", sum(all_results$qc_pass), "/", nrow(all_results))
message("IVW results written: ", file.path(result_dir, "finngen_primary_ivw_summary.tsv"))
message("Error log written: ", error_file, " (", length(errors), " entries)")
if (any(!all_results$qc_pass)) quit(status = 1L)
