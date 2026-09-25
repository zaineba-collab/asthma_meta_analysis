#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

# Standard sensitivity analyses only: Asthma -> 12 local FinnGen outcomes.
# This script does not harmonise, clump, filter instruments, or run additional MR methods.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(arg)) sub("^--file=", "", arg[1]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)

ready_dir <- p("results", "mr", "finngen", "harmonised", "mr_ready")
config_file <- p("config", "finngen_mr_outcomes.tsv")
harm_qc_file <- p("results", "mr", "finngen", "qc", "finngen_harmonisation_qc.tsv")
primary_file <- p("results", "mr", "finngen", "results", "primary", "finngen_primary_mr_all_methods.tsv")
sens_dir <- p("results", "mr", "finngen", "sensitivity")
loo_dir <- file.path(sens_dir, "leave_one_out")
single_dir <- file.path(sens_dir, "single_snp")
plot_root <- p("results", "mr", "finngen", "plots")
funnel_dir <- file.path(plot_root, "funnel")
loo_plot_dir <- file.path(plot_root, "leave_one_out")
scatter_dir <- file.path(plot_root, "scatter")
forest_dir <- file.path(plot_root, "forest")
qc_dir <- p("results", "mr", "finngen", "qc")
error_file <- file.path(qc_dir, "finngen_sensitivity_errors.tsv")
for (d in c(sens_dir, loo_dir, single_dir, funnel_dir, loo_plot_dir, scatter_dir, forest_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

for (f in c(config_file, harm_qc_file, primary_file)) if (!file.exists(f)) stop("Missing required input: ", f, call. = FALSE)
cfg <- fread(config_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
harm_qc <- fread(harm_qc_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
primary <- fread(primary_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
if (nrow(cfg) != 12L || uniqueN(cfg$outcome_id) != 12L) stop("Outcome config must contain 12 unique outcomes", call. = FALSE)
if (nrow(primary) != 60L) stop("Primary MR table must contain 60 estimates", call. = FALSE)

required <- c("SNP", "mr_keep", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome",
              "effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome",
              "other_allele.outcome", "exposure", "outcome", "id.exposure", "id.outcome")
errors <- list(); heterogeneity <- list(); egger <- list(); loo_summary <- list(); single_summary <- list()
add_error <- function(outcome, outcome_id, diagnostic, message) {
  errors[[length(errors) + 1L]] <<- data.table(
    outcome = as.character(outcome), outcome_id = as.character(outcome_id),
    diagnostic = diagnostic, error = as.character(message)
  )
}
save_plot <- function(plot, pdf_file, png_file, outcome, outcome_id, diagnostic, width = 8, height = 6) {
  tryCatch({
    ggsave(pdf_file, plot = plot, width = width, height = height, units = "in", limitsize = FALSE)
    if (diagnostic %in% c("scatter_plot", "funnel_plot", "leave_one_out_plot")) {
      plot <- plot + labs(title = as.character(outcome)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    }
    ggsave(png_file, plot = plot, width = width, height = height, units = "in", dpi = 300, limitsize = FALSE)
  }, error = function(e) add_error(outcome, outcome_id, diagnostic, conditionMessage(e)))
}
finite_probability <- function(x) all(is.finite(x)) && all(x >= 0 & x <= 1)

for (i in seq_len(nrow(cfg))) {
  meta <- cfg[i]
  input_file <- file.path(ready_dir, paste0(meta$output_stub, "_mr_ready.tsv"))
  message("Sensitivity analyses: Asthma -> ", meta$outcome)
  if (!file.exists(input_file)) {
    add_error(meta$outcome, meta$outcome_id, "input_validation", paste("Missing MR-ready input:", input_file))
    next
  }
  dat <- tryCatch(fread(input_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id")),
                  error = function(e) { add_error(meta$outcome, meta$outcome_id, "input_validation", conditionMessage(e)); NULL })
  if (is.null(dat)) next
  validation_error <- tryCatch({
    absent <- setdiff(required, names(dat)); if (length(absent)) stop("Missing fields: ", paste(absent, collapse = ", "))
    if (anyDuplicated(dat$SNP)) stop("SNP is not unique")
    if (!nrow(dat) || anyNA(dat$mr_keep) || any(!dat$mr_keep)) stop("Not all rows have mr_keep == TRUE")
    fields <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
    if (!all(vapply(dat[, ..fields], is.numeric, logical(1)))) stop("Beta/SE fields are not numeric")
    if (anyNA(dat[, ..fields]) || any(!is.finite(as.matrix(dat[, ..fields])))) stop("Missing or non-finite beta/SE")
    if (any(dat$se.exposure <= 0) || any(dat$se.outcome <= 0)) stop("SE is non-positive")
    if (anyNA(dat[, .(effect_allele.exposure, other_allele.exposure, effect_allele.outcome, other_allele.outcome)]) ||
        any(dat$effect_allele.exposure != dat$effect_allele.outcome) ||
        any(dat$other_allele.exposure != dat$other_allele.outcome)) stop("Alleles are not harmonised")
    expected <- harm_qc[outcome_id == as.character(meta$outcome_id), final_mr_keep_snps]
    if (length(expected) != 1L || nrow(dat) != expected) stop("SNP count does not match harmonisation QC")
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(validation_error)) {
    add_error(meta$outcome, meta$outcome_id, "input_validation", validation_error)
    next
  }
  n_instruments <- nrow(dat)
  ivw <- primary[outcome_id == as.character(meta$outcome_id) & method_id == "mr_ivw"]
  if (nrow(ivw) != 1L || !is.finite(ivw$beta) || !is.finite(ivw$p_value)) {
    add_error(meta$outcome, meta$outcome_id, "primary_result", "Missing or invalid saved IVW estimate")
    next
  }

  het <- tryCatch(as.data.table(TwoSampleMR::mr_heterogeneity(as.data.frame(dat),
    method_list = c("mr_ivw", "mr_egger_regression"))),
    error = function(e) { add_error(meta$outcome, meta$outcome_id, "heterogeneity", conditionMessage(e)); NULL })
  if (!is.null(het)) {
    het[, `:=`(outcome_id = as.character(meta$outcome_id), outcome_type = as.character(meta$outcome_type), nsnp = n_instruments)]
    het[, method_id := fifelse(method == "Inverse variance weighted", "mr_ivw",
                               fifelse(method == "MR Egger", "mr_egger_regression", NA_character_))]
    het <- het[, .(exposure, outcome, outcome_id, outcome_type, method, method_id, nsnp,
                   Q, Q_df, Q_p_value = Q_pval)]
    if (nrow(het) != 2L || anyNA(het$method_id) || any(!is.finite(het$Q)) ||
        any(het$Q_df < 0) || !finite_probability(het$Q_p_value)) {
      add_error(meta$outcome, meta$outcome_id, "heterogeneity_qc", "Invalid or incomplete heterogeneity results")
    }
    fwrite(het, file.path(sens_dir, paste0(meta$output_stub, "_heterogeneity.tsv")), sep = "\t", quote = FALSE, na = "NA")
    heterogeneity[[i]] <- het
  }

  eg <- tryCatch(as.data.table(TwoSampleMR::mr_pleiotropy_test(as.data.frame(dat))),
    error = function(e) { add_error(meta$outcome, meta$outcome_id, "egger_intercept", conditionMessage(e)); NULL })
  if (!is.null(eg)) {
    eg <- eg[, .(exposure, outcome, outcome_id = as.character(meta$outcome_id),
                 outcome_type = as.character(meta$outcome_type), nsnp = n_instruments,
                 egger_intercept, SE = se, p_value = pval)]
    if (nrow(eg) != 1L || any(!is.finite(as.matrix(eg[, .(egger_intercept, SE)]))) || eg$SE <= 0 || !finite_probability(eg$p_value)) {
      add_error(meta$outcome, meta$outcome_id, "egger_intercept_qc", "Invalid Egger-intercept result")
    }
    egger[[i]] <- eg
  }

  loo <- tryCatch(as.data.table(TwoSampleMR::mr_leaveoneout(as.data.frame(dat))),
    error = function(e) { add_error(meta$outcome, meta$outcome_id, "leave_one_out", conditionMessage(e)); NULL })
  if (!is.null(loo)) {
    setnames(loo, c("b", "se", "p"), c("beta", "SE", "p_value"))
    loo[, `:=`(outcome_id = as.character(meta$outcome_id), outcome_type = as.character(meta$outcome_type))]
    setcolorder(loo, c("exposure", "outcome", "outcome_id", "outcome_type", "SNP", "beta", "SE", "p_value",
                       setdiff(names(loo), c("exposure", "outcome", "outcome_id", "outcome_type", "SNP", "beta", "SE", "p_value"))))
    excluded <- loo[SNP != "All"]
    if (nrow(excluded) != n_instruments || anyDuplicated(excluded$SNP) ||
        any(!dat$SNP %in% excluded$SNP) || anyNA(excluded[, .(beta, SE, p_value)])) {
      add_error(meta$outcome, meta$outcome_id, "leave_one_out_qc", "Leave-one-out rows do not represent every unique instrument")
    }
    changes <- abs(excluded$beta - ivw$beta)
    influential_index <- which.max(changes)
    sign_change <- any(sign(excluded$beta) != sign(ivw$beta))
    loo_summary[[i]] <- data.table(
      outcome = as.character(meta$outcome), outcome_id = as.character(meta$outcome_id), nsnp = n_instruments,
      original_ivw_beta = ivw$beta, minimum_leave_one_out_beta = min(excluded$beta),
      maximum_leave_one_out_beta = max(excluded$beta), minimum_leave_one_out_p_value = min(excluded$p_value),
      maximum_leave_one_out_p_value = max(excluded$p_value), leave_one_out_sign_change = sign_change,
      most_influential_snp = excluded$SNP[influential_index],
      maximum_absolute_beta_change = changes[influential_index]
    )
    fwrite(loo, file.path(loo_dir, paste0(meta$output_stub, "_leave_one_out.tsv")), sep = "\t", quote = FALSE, na = "NA")
    loo_plot <- tryCatch(TwoSampleMR::mr_leaveoneout_plot(as.data.frame(loo[, .(
      exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b = beta, se = SE, p = p_value)]))[[1]],
      error = function(e) { add_error(meta$outcome, meta$outcome_id, "leave_one_out_plot", conditionMessage(e)); NULL })
    if (!is.null(loo_plot)) save_plot(loo_plot,
      file.path(loo_plot_dir, paste0(meta$output_stub, "_leave_one_out.pdf")),
      file.path(loo_plot_dir, paste0(meta$output_stub, "_leave_one_out.png")),
      meta$outcome, meta$outcome_id, "leave_one_out_plot")
  }

  single <- tryCatch(as.data.table(TwoSampleMR::mr_singlesnp(as.data.frame(dat))),
    error = function(e) { add_error(meta$outcome, meta$outcome_id, "single_snp", conditionMessage(e)); NULL })
  if (!is.null(single)) {
    setnames(single, c("b", "se", "p"), c("beta", "SE", "p_value"))
    single[, `:=`(
      outcome_id = as.character(meta$outcome_id), outcome_type = as.character(meta$outcome_type),
      method = fifelse(grepl("^All - ", SNP), sub("^All - ", "", SNP), "Wald ratio")
    )]
    setcolorder(single, c("exposure", "outcome", "outcome_id", "outcome_type", "SNP", "method", "beta", "SE", "p_value",
                          setdiff(names(single), c("exposure", "outcome", "outcome_id", "outcome_type", "SNP", "method", "beta", "SE", "p_value"))))
    individual <- single[!grepl("^All - ", SNP)]
    if (nrow(individual) != n_instruments || anyDuplicated(individual$SNP) ||
        any(!dat$SNP %in% individual$SNP) || anyNA(individual[, .(beta, SE, p_value)])) {
      add_error(meta$outcome, meta$outcome_id, "single_snp_qc", "Single-SNP rows do not represent every unique instrument")
    }
    strongest_p <- which.min(individual$p_value)
    largest_ratio <- which.max(abs(individual$beta))
    single_summary[[i]] <- data.table(
      outcome = as.character(meta$outcome), outcome_id = as.character(meta$outcome_id), nsnp = n_instruments,
      strongest_single_snp = individual$SNP[strongest_p], strongest_single_snp_p_value = individual$p_value[strongest_p],
      largest_absolute_wald_ratio_snp = individual$SNP[largest_ratio],
      largest_absolute_wald_ratio = abs(individual$beta[largest_ratio])
    )
    fwrite(single, file.path(single_dir, paste0(meta$output_stub, "_single_snp.tsv")), sep = "\t", quote = FALSE, na = "NA")

    plot_single <- single[, .(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b = beta, se = SE, p = p_value)]
    funnel_plot <- tryCatch(TwoSampleMR::mr_funnel_plot(as.data.frame(plot_single))[[1]],
      error = function(e) { add_error(meta$outcome, meta$outcome_id, "funnel_plot", conditionMessage(e)); NULL })
    if (!is.null(funnel_plot)) save_plot(funnel_plot,
      file.path(funnel_dir, paste0(meta$output_stub, "_funnel.pdf")), file.path(funnel_dir, paste0(meta$output_stub, "_funnel.png")),
      meta$outcome, meta$outcome_id, "funnel_plot")

    # Complete numerical results are retained above. A readable forest plot uses
    # the 30 strongest individual associations plus the aggregate IVW/Egger rows.
    forest_subset <- rbind(single[!grepl("^All - ", SNP)][order(p_value)][seq_len(min(30L, .N))],
                           single[grepl("^All - ", SNP)], fill = TRUE)
    forest_input <- forest_subset[, .(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b = beta, se = SE, p = p_value)]
    forest_plot <- tryCatch(TwoSampleMR::mr_forest_plot(as.data.frame(forest_input))[[1]] +
      labs(title = paste0(meta$outcome, ": 30 strongest SNPs (subset)")),
      error = function(e) { add_error(meta$outcome, meta$outcome_id, "forest_plot", conditionMessage(e)); NULL })
    if (!is.null(forest_plot)) save_plot(forest_plot,
      file.path(forest_dir, paste0(meta$output_stub, "_forest_top30_subset.pdf")),
      file.path(forest_dir, paste0(meta$output_stub, "_forest_top30_subset.png")),
      meta$outcome, meta$outcome_id, "forest_plot", width = 8, height = 10)
  }

  # Reformat the saved primary estimates into the native result schema required
  # by mr_scatter_plot; no primary estimator is rerun here.
  saved <- primary[outcome_id == as.character(meta$outcome_id)]
  scatter_result <- saved[, .(id.exposure = "asthma_gwama", id.outcome = as.character(meta$outcome_id),
    outcome, exposure, method, nsnp, b = beta, se = SE, pval = p_value)]
  scatter_plot <- tryCatch(TwoSampleMR::mr_scatter_plot(as.data.frame(scatter_result), as.data.frame(dat))[[1]],
    error = function(e) { add_error(meta$outcome, meta$outcome_id, "scatter_plot", conditionMessage(e)); NULL })
  if (!is.null(scatter_plot)) save_plot(scatter_plot,
    file.path(scatter_dir, paste0(meta$output_stub, "_scatter.pdf")), file.path(scatter_dir, paste0(meta$output_stub, "_scatter.png")),
    meta$outcome, meta$outcome_id, "scatter_plot")
}

het_all <- rbindlist(heterogeneity, fill = TRUE)
egger_all <- rbindlist(egger, fill = TRUE)
loo_all <- rbindlist(loo_summary, fill = TRUE)
single_all <- rbindlist(single_summary, fill = TRUE)
for (x in list(het_all, egger_all, loo_all, single_all)) {
  if (nrow(x)) x[, outcome_order := match(outcome_id, cfg$outcome_id)]
}
if (nrow(het_all)) setorder(het_all, outcome_order, method_id)
if (nrow(egger_all)) setorder(egger_all, outcome_order)
if (nrow(loo_all)) setorder(loo_all, outcome_order)
if (nrow(single_all)) setorder(single_all, outcome_order)
egger_all[, p_egger_fdr := p.adjust(p_value, method = "BH")]
egger_all[, directional_pleiotropy_evidence := p_value < 0.05]
fwrite(het_all[, !"outcome_order"], file.path(sens_dir, "finngen_mr_heterogeneity.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(egger_all[, !"outcome_order"], file.path(sens_dir, "finngen_mr_egger_intercept.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(loo_all[, !"outcome_order"], file.path(sens_dir, "finngen_leave_one_out_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(single_all[, !"outcome_order"], file.path(sens_dir, "finngen_single_snp_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

ivw_primary <- primary[method_id == "mr_ivw", .(outcome, outcome_id, nsnp, ivw_beta = beta, ivw_p_value = p_value)]
ivw_het <- het_all[method_id == "mr_ivw", .(outcome_id, ivw_Q = Q, ivw_Q_df = Q_df, ivw_heterogeneity_p_value = Q_p_value)]
egger_het <- het_all[method_id == "mr_egger_regression", .(outcome_id, mr_egger_Q = Q, mr_egger_Q_df = Q_df,
                                                            mr_egger_heterogeneity_p_value = Q_p_value)]
summary <- Reduce(function(x, y) merge(x, y, by = "outcome_id", all = TRUE, sort = FALSE), list(
  ivw_primary, ivw_het, egger_het,
  egger_all[, .(outcome_id, egger_intercept, egger_intercept_SE = SE, egger_intercept_p_value = p_value,
                p_egger_fdr, directional_pleiotropy_evidence)],
  loo_all[, .(outcome_id, leave_one_out_sign_change, most_influential_leave_one_out_snp = most_influential_snp,
              maximum_absolute_leave_one_out_beta_change = maximum_absolute_beta_change)],
  single_all[, .(outcome_id, strongest_single_snp, strongest_single_snp_p_value)]
))
summary[, heterogeneity_evidence := ivw_heterogeneity_p_value < 0.05]
summary[, outcome_order := match(outcome_id, cfg$outcome_id)]
setorder(summary, outcome_order)
setcolorder(summary, c("outcome", "outcome_id", "nsnp", setdiff(names(summary), c("outcome", "outcome_id", "nsnp", "outcome_order")), "outcome_order"))
fwrite(summary[, !"outcome_order"], file.path(sens_dir, "finngen_mr_sensitivity_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

expected_plot_files <- unlist(lapply(c(funnel_dir, loo_plot_dir, scatter_dir, forest_dir), function(d) {
  c(file.path(d, paste0(cfg$output_stub, if (d == funnel_dir) "_funnel.pdf" else if (d == loo_plot_dir) "_leave_one_out.pdf" else if (d == scatter_dir) "_scatter.pdf" else "_forest_top30_subset.pdf")),
    file.path(d, paste0(cfg$output_stub, if (d == funnel_dir) "_funnel.png" else if (d == loo_plot_dir) "_leave_one_out.png" else if (d == scatter_dir) "_scatter.png" else "_forest_top30_subset.png")))
}))
missing_plots <- expected_plot_files[!file.exists(expected_plot_files) | file.info(expected_plot_files)$size == 0]
if (length(missing_plots)) for (f in missing_plots) add_error(NA, NA, "plot_validation", paste("Missing or empty plot:", f))
if (nrow(summary) != 12L) add_error(NA, NA, "combined_qc", paste("Sensitivity summary has", nrow(summary), "rather than 12 outcomes"))

if (length(errors)) {
  fwrite(rbindlist(errors, fill = TRUE), error_file, sep = "\t", quote = FALSE, na = "NA")
} else {
  fwrite(data.table(outcome = character(), outcome_id = character(), diagnostic = character(), error = character()),
         error_file, sep = "\t", quote = FALSE, na = "NA")
}
message("Outcomes in sensitivity summary: ", nrow(summary), "/12")
message("Plots generated: ", sum(file.exists(expected_plot_files)), "/", length(expected_plot_files))
message("Error log: ", error_file, " (", length(errors), " entries)")
if (length(errors)) quit(status = 1L)
