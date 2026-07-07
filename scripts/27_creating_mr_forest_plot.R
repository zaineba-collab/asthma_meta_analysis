#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)

# Read the reproducible IVW summary created by script 26.
mr_plot <- fread("results/mr/summary/asthma_multi_outcome_ivw_summary.tsv") %>%
  as.data.frame() %>%
  mutate(
    outcome = outcome_name,
    beta = beta_estimate,
    se = standard_error,
    pval = p_value,
    lower = beta_lci95,
    upper = beta_uci95,
    outcome = factor(outcome, levels = rev(outcome))
  )

# Create a publication-style forest plot of the main IVW estimates.
forest_plot <- ggplot(mr_plot, aes(x = beta, y = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  labs(
    x = "MR effect estimate",
    y = NULL,
    title = "Mendelian randomisation estimates for asthma and cardiometabolic traits",
    subtitle = "Main IVW estimates with 95% confidence intervals"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 15),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 17)
  )

ggsave(
  "results/mr/summary/asthma_multi_outcome_ivw_forest_plot_classic.png",
  plot = forest_plot,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  "results/mr/summary/asthma_multi_outcome_ivw_forest_plot_classic.pdf",
  plot = forest_plot,
  width = 10,
  height = 6
)

message("Written classic MR forest plot PNG and PDF.")
