# Copyright 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# https://www.apache.org/licenses/LICENSE-2.0
#
# The code is provided "AS IS", without guarantees or warranties of any kind.
# To the extent permitted by law, the copyright holder accepts no responsibility or liability for use of the code or its outputs.
# See LICENSE for the full warranty disclaimer and limitation of liability.
# Run from the repository root; generated files are written in the working directory.
# Minimal GO dot plot from PANTHER Excel -----------------------

# install.packages(c("readxl","dplyr","ggplot2","stringr","janitor","forcats","scales"))
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(janitor)
library(forcats)
library(scales)

################################################################################
# DESeq2 mRNA 6m
################################################################################

# 1) Read Excel
file  <- "PANTHER_GO-Slim_DESeq2_mRNA_6m.xlsx"   # <-- your PANTHER GO-Slim BP results
sheet <- 1
df_raw <- read_excel(file, sheet = sheet)
df <- janitor::clean_names(df_raw)

# 2) Build tidy table
observed_vec <- df %>%
  dplyr::select(starts_with("dysregulated_genes_")) %>%  # first dysregulated_genes_* col
  dplyr::pull(1)

plot_df <- tibble(
  term_full    = df$panther_go_slim_biological_process,
  go_id        = stringr::str_extract(df$panther_go_slim_biological_process, "GO:\\d{7}"),
  term         = stringr::str_trim(stringr::str_remove(df$panther_go_slim_biological_process, "\\s*\\(GO:\\d{7}\\)")),
  observed     = observed_vec,
  expected     = as.numeric(df$dysregulated_genes_expected),         # force numeric
  fold_enrich  = as.numeric(df$dysregulated_genes_fold_enrichment),  # force numeric
  fdr          = suppressWarnings(as.numeric(df$dysregulated_genes_fdr))
) %>%
  dplyr::filter(!is.na(fdr), is.finite(fdr), fdr > 0) %>%
  dplyr::mutate(neg_log10_fdr = -log10(fdr))

# 3) Pick how many terms to show
top_n <- 30
df_top <- plot_df %>%
  dplyr::filter(fold_enrich > 1) %>%                         # shows only overrepresented terms, not underrepresented
  dplyr::arrange(fdr, dplyr::desc(fold_enrich)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::mutate(
    term_label = stringr::str_wrap(term, width = 50),
    term_label = forcats::fct_reorder(term_label, neg_log10_fdr)  # order by significance
  )

# 4) Dot plot
p_mRNA_6m <- ggplot(df_top, aes(x = fold_enrich, y = term_label)) +
  geom_point(aes(size = observed, color = neg_log10_fdr)) +
  scale_size_area(max_size = 10, name = "Observed genes") +
  scale_x_continuous(
    limits = c(1, NA),  # force axis to start at 1; let ggplot set the upper limit
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 6)  # ~6 evenly spaced ticks
  ) +
  labs(
    x = "Fold enrichment",
    y = NULL,
    color = expression(-log[10]~FDR),
    title = "Overrepresented GO terms at 6 months",
    subtitle = paste0("Top ", top_n, " terms by FDR")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

print(p_mRNA_6m)

# 5) Save
ggsave("GO_dotplot_TDP_mRNA_6m.svg", plot = p_mRNA_6m, width = 8, height = 7, device = "svg")

################################################################################
# DESeq2 mRNA 9m
################################################################################

# 1) Read Excel
file  <- "PANTHER_GO-Slim_DESeq2_mRNA_9m.xlsx"   # <-- your PANTHER GO-Slim BP results
sheet <- 1
df_raw <- read_excel(file, sheet = sheet)
df <- janitor::clean_names(df_raw)

# 2) Build tidy table
observed_vec <- df %>%
  dplyr::select(starts_with("dysregulated_genes_")) %>%  # first dysregulated_genes_* col
  dplyr::pull(1)

plot_df <- tibble(
  term_full    = df$panther_go_slim_biological_process,
  go_id        = stringr::str_extract(df$panther_go_slim_biological_process, "GO:\\d{7}"),
  term         = stringr::str_trim(stringr::str_remove(df$panther_go_slim_biological_process, "\\s*\\(GO:\\d{7}\\)")),
  observed     = observed_vec,
  expected     = as.numeric(df$dysregulated_genes_expected),         # force numeric
  fold_enrich  = as.numeric(df$dysregulated_genes_fold_enrichment),  # force numeric
  fdr          = suppressWarnings(as.numeric(df$dysregulated_genes_fdr))
) %>%
  dplyr::filter(!is.na(fdr), is.finite(fdr), fdr > 0) %>%
  dplyr::mutate(neg_log10_fdr = -log10(fdr))

# 3) Pick how many terms to show
top_n <- 30
df_top <- plot_df %>%
  dplyr::filter(fold_enrich > 1) %>%                         # shows only overrepresented terms, not underrepresented
  dplyr::arrange(fdr, dplyr::desc(fold_enrich)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::mutate(
    term_label = stringr::str_wrap(term, width = 50),
    term_label = forcats::fct_reorder(term_label, neg_log10_fdr)  # order by significance
  )

# 4) Dot plot
p_mRNA_9m <- ggplot(df_top, aes(x = fold_enrich, y = term_label)) +
  geom_point(aes(size = observed, color = neg_log10_fdr)) +
  scale_size_area(max_size = 10, name = "Observed genes") +
  scale_x_continuous(
    limits = c(1, NA),  # force axis to start at 1; let ggplot set the upper limit
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 6)  # ~6 evenly spaced ticks
  ) +
  labs(
    x = "Fold enrichment",
    y = NULL,
    color = expression(-log[10]~FDR),
    title = "Overrepresented GO terms at 9 months",
    subtitle = paste0("Top ", top_n, " terms by FDR")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

print(p_mRNA_9m)

# 5) Save
ggsave("GO_dotplot_TDP_mRNA_9m.svg", plot = p_mRNA_9m, width = 8, height = 7, device = "svg")

################################################################################
# miRPath miRNA targets 6m
################################################################################

# 1) Read Excel
file  <- "PANTHER_GO-Slim_miRPath_miRNA_6m.xlsx"
sheet <- 1
df_raw <- read_excel(file, sheet = sheet)
df <- janitor::clean_names(df_raw)

# Helper: safely pull a column by any of several candidate names
get_col <- function(d, candidates) {
  hit <- intersect(candidates, names(d))
  if (length(hit)) d[[hit[1]]] else NULL
}

# Derive mirna count if needed from "miRNA Names"
derive_miRNA_n <- function(names_col, n_rows) {
  if (is.null(names_col)) return(rep(NA_integer_, n_rows))
  sapply(stringr::str_split(as.character(names_col), "[;,]"), function(x) {
    sum(nzchar(trimws(x)))
  })
}

# 2) Build tidy table (robust)
term_full   <- get_col(df, c("go_biological_process_term_name", "term_name", "go_term"))
term_genes  <- suppressWarnings(as.numeric(get_col(df, c("term_genes", "genes_in_term", "term_genes_n"))))
targets_n   <- suppressWarnings(as.numeric(get_col(df, c("target_genes_n", "targets_n", "target_genes"))))
mirnas_n    <- suppressWarnings(as.numeric(get_col(df, c("mirnas_n", "mirna_n", "n_mirnas"))))
mirna_names <- get_col(df, c("mirna_names", "mirnas", "mirna_list"))
p_sep       <- suppressWarnings(as.numeric(get_col(df, c("p_value_separate", "p_value", "p_separate"))))
merged_p    <- suppressWarnings(as.numeric(get_col(df, c("merged_p_value", "p_merged", "combined_p_value"))))
merged_fdr  <- suppressWarnings(as.numeric(get_col(df, c("merged_fdr", "fdr_merged", "adjusted_p_value", "adj_p"))))

# If miRNA count missing or length-0, derive from names
if (is.null(mirnas_n) || length(mirnas_n) == 0) {
  mirnas_n <- derive_miRNA_n(mirna_names, nrow(df))
}


# --- after you've read and cleaned df_raw -> df ---

# Parse "Target Genes (n)" which may be a comma/semicolon-separated list like "3,3,4"
targets_col_raw <- get_col(df, c("target_genes_n", "targets_n", "target_genes"))

parse_counts <- function(x) {
  if (is.null(x)) return(rep(NA_real_, nrow(df)))
  sapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    vals <- suppressWarnings(as.numeric(unlist(stringr::str_split(s, "[;,]"))))
    # choose an aggregator: max is conservative for "any miRNA hits this many genes"
    if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
  })
}

targets_n_agg <- parse_counts(targets_col_raw)

# Keep your other columns as before...
term_full   <- get_col(df, c("go_biological_process_term_name", "term_name", "go_term"))
term_genes  <- suppressWarnings(as.numeric(get_col(df, c("term_genes", "genes_in_term", "term_genes_n"))))
mirnas_n    <- suppressWarnings(as.numeric(get_col(df, c("mirnas_n", "mirna_n", "n_mirnas"))))
mirna_names <- get_col(df, c("mirna_names", "mirnas", "mirna_list"))
merged_p    <- suppressWarnings(as.numeric(get_col(df, c("merged_p_value", "p_merged", "combined_p_value"))))
merged_fdr  <- suppressWarnings(as.numeric(get_col(df, c("merged_fdr", "fdr_merged", "adjusted_p_value", "adj_p"))))

# Derive miRNA count if missing
if (is.null(mirnas_n) || length(mirnas_n) == 0) {
  mirnas_n <- derive_miRNA_n(mirna_names, nrow(df))
}

plot_df <- tibble::tibble(
  term_full   = as.character(term_full),
  term        = stringr::str_trim(as.character(term_full)),
  term_genes  = term_genes,
  targets_n   = targets_n_agg,                # use parsed/aggregated counts
  mirnas_n    = mirnas_n,
  merged_p    = merged_p,
  merged_fdr  = merged_fdr
) |>
  dplyr::filter(!is.na(merged_fdr), is.finite(merged_fdr), merged_fdr > 0) |>
  dplyr::mutate(
    neg_log10_fdr = -log10(merged_fdr),
    hit_fraction  = ifelse(term_genes > 0, targets_n / term_genes, NA_real_)
  )

# Keep only plottable rows
plotable_df <- plot_df |>
  dplyr::filter(
    !is.na(term_full),
    is.finite(neg_log10_fdr),
    !is.na(targets_n), is.finite(targets_n)
  )

# Choose top 30 by smallest FDR, then by larger targets_n
top_n <- 30
df_top <- plotable_df %>%
  arrange(merged_fdr, desc(targets_n)) %>%   # 1) sort by FDR then targets
  slice_head(n = top_n) %>%
  mutate(
    term_label = stringr::str_wrap(term_full, width = 50),
    # 2) reorder factor explicitly by merged_fdr
    term_label = forcats::fct_reorder(term_label, merged_fdr, .desc = TRUE)
  )


# Primary plot: x = -log10(FDR)  (recommended)
p_sig <- ggplot2::ggplot(df_top, aes(x = neg_log10_fdr, y = term_label)) +
  ggplot2::geom_point(aes(size = targets_n, color = neg_log10_fdr)) +
  ggplot2::scale_size_area(max_size = 10, limits = c(0, NA), name = "Observed genes") +
  ggplot2::scale_x_continuous(
    limits = c(1, NA),
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  ggplot2::labs(
    x = expression(-log[10]~"(p-value)"),
    y = NULL,
    color = expression(-log[10]~"(p-value)"),
    title = "Overrepresented GO terms at 6 months",
    subtitle = paste0("Top ", top_n, " terms by unadjusted p-value")
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title  = ggplot2::element_text(face = "bold")
  )

print(p_sig)
ggplot2::ggsave("GO_dotplot_TDP_miRNA_6m.svg", plot = p_sig, width = 8, height = 9, device = "svg")

################################################################################
# miRPath miRNA targets 9m
################################################################################

# 1) Read Excel
file  <- "PANTHER_GO-Slim_miRPath_miRNA_9m.xlsx"
sheet <- 1
df_raw <- read_excel(file, sheet = sheet)
df <- janitor::clean_names(df_raw)

# Helper: safely pull a column by any of several candidate names
get_col <- function(d, candidates) {
  hit <- intersect(candidates, names(d))
  if (length(hit)) d[[hit[1]]] else NULL
}

# Derive mirna count if needed from "miRNA Names"
derive_miRNA_n <- function(names_col, n_rows) {
  if (is.null(names_col)) return(rep(NA_integer_, n_rows))
  sapply(stringr::str_split(as.character(names_col), "[;,]"), function(x) {
    sum(nzchar(trimws(x)))
  })
}

# 2) Build tidy table (robust)
term_full   <- get_col(df, c("go_biological_process_term_name", "term_name", "go_term"))
term_genes  <- suppressWarnings(as.numeric(get_col(df, c("term_genes", "genes_in_term", "term_genes_n"))))
targets_n   <- suppressWarnings(as.numeric(get_col(df, c("target_genes_n", "targets_n", "target_genes"))))
mirnas_n    <- suppressWarnings(as.numeric(get_col(df, c("mirnas_n", "mirna_n", "n_mirnas"))))
mirna_names <- get_col(df, c("mirna_names", "mirnas", "mirna_list"))
p_sep       <- suppressWarnings(as.numeric(get_col(df, c("p_value_separate", "p_value", "p_separate"))))
merged_p    <- suppressWarnings(as.numeric(get_col(df, c("merged_p_value", "p_merged", "combined_p_value"))))
merged_fdr  <- suppressWarnings(as.numeric(get_col(df, c("merged_fdr", "fdr_merged", "adjusted_p_value", "adj_p"))))

# If miRNA count missing or length-0, derive from names
if (is.null(mirnas_n) || length(mirnas_n) == 0) {
  mirnas_n <- derive_miRNA_n(mirna_names, nrow(df))
}


# --- after you've read and cleaned df_raw -> df ---

# Parse "Target Genes (n)" which may be a comma/semicolon-separated list like "3,3,4"
targets_col_raw <- get_col(df, c("target_genes_n", "targets_n", "target_genes"))

parse_counts <- function(x) {
  if (is.null(x)) return(rep(NA_real_, nrow(df)))
  sapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    vals <- suppressWarnings(as.numeric(unlist(stringr::str_split(s, "[;,]"))))
    # choose an aggregator: max is conservative for "any miRNA hits this many genes"
    if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
  })
}

targets_n_agg <- parse_counts(targets_col_raw)

# Keep your other columns as before...
term_full   <- get_col(df, c("go_biological_process_term_name", "term_name", "go_term"))
term_genes  <- suppressWarnings(as.numeric(get_col(df, c("term_genes", "genes_in_term", "term_genes_n"))))
mirnas_n    <- suppressWarnings(as.numeric(get_col(df, c("mirnas_n", "mirna_n", "n_mirnas"))))
mirna_names <- get_col(df, c("mirna_names", "mirnas", "mirna_list"))
merged_p    <- suppressWarnings(as.numeric(get_col(df, c("merged_p_value", "p_merged", "combined_p_value"))))
merged_fdr  <- suppressWarnings(as.numeric(get_col(df, c("merged_fdr", "fdr_merged", "adjusted_p_value", "adj_p"))))

# Derive miRNA count if missing
if (is.null(mirnas_n) || length(mirnas_n) == 0) {
  mirnas_n <- derive_miRNA_n(mirna_names, nrow(df))
}

plot_df <- tibble::tibble(
  term_full   = as.character(term_full),
  term        = stringr::str_trim(as.character(term_full)),
  term_genes  = term_genes,
  targets_n   = targets_n_agg,                # use parsed/aggregated counts
  mirnas_n    = mirnas_n,
  merged_p    = merged_p,
  merged_fdr  = merged_fdr
) |>
  dplyr::filter(!is.na(merged_fdr), is.finite(merged_fdr), merged_fdr > 0) |>
  dplyr::mutate(
    neg_log10_fdr = -log10(merged_fdr),
    hit_fraction  = ifelse(term_genes > 0, targets_n / term_genes, NA_real_)
  )

# Keep only plottable rows
plotable_df <- plot_df |>
  dplyr::filter(
    !is.na(term_full),
    is.finite(neg_log10_fdr),
    !is.na(targets_n), is.finite(targets_n)
  )

# Choose top 30 by smallest FDR, then by larger targets_n
top_n <- 30
df_top <- plotable_df %>%
  arrange(merged_fdr, desc(targets_n)) %>%   # 1) sort by FDR then targets
  slice_head(n = top_n) %>%
  mutate(
    term_label = stringr::str_wrap(term_full, width = 50),
    # 2) reorder factor explicitly by merged_fdr
    term_label = forcats::fct_reorder(term_label, merged_fdr, .desc = TRUE)
  )


# Primary plot: x = -log10(FDR)  (recommended)
p_sig <- ggplot2::ggplot(df_top, aes(x = neg_log10_fdr, y = term_label)) +
  ggplot2::geom_point(aes(size = targets_n, color = neg_log10_fdr)) +
  ggplot2::scale_size_area(max_size = 10, limits = c(0, NA), name = "Observed genes") +
  ggplot2::scale_x_continuous(
    limits = c(1, NA),
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  ggplot2::labs(
    x = expression(-log[10]~"(p-value)"),
    y = NULL,
    color = expression(-log[10]~"(p-value)"),
    title = "Overrepresented GO terms at 9 months",
    subtitle = paste0("Top ", top_n, " terms by unadjusted p-value")
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title  = ggplot2::element_text(face = "bold")
  )

print(p_sig)
ggplot2::ggsave("GO_dotplot_TDP_miRNA_9m.svg", plot = p_sig, width = 8, height = 5, device = "svg")


################################################################################
# miRPath miRNA targets 12m
################################################################################

# 1) Read Excel
file  <- "PANTHER_GO-Slim_miRPath_miRNA_12m.xlsx"
sheet <- 1
df_raw <- read_excel(file, sheet = sheet)
df <- janitor::clean_names(df_raw)

# Helper: safely pull a column by any of several candidate names
get_col <- function(d, candidates) {
  hit <- intersect(candidates, names(d))
  if (length(hit)) d[[hit[1]]] else NULL
}

# Derive mirna count if needed from "miRNA Names"
derive_miRNA_n <- function(names_col, n_rows) {
  if (is.null(names_col)) return(rep(NA_integer_, n_rows))
  sapply(stringr::str_split(as.character(names_col), "[;,]"), function(x) {
    sum(nzchar(trimws(x)))
  })
}

# 2) Build tidy table (robust)
term_full   <- get_col(df, c("go_biological_process_term_name", "term_name", "go_term"))
term_genes  <- suppressWarnings(as.numeric(get_col(df, c("term_genes", "genes_in_term", "term_genes_n"))))
targets_n   <- suppressWarnings(as.numeric(get_col(df, c("target_genes_n", "targets_n", "target_genes"))))
mirnas_n    <- suppressWarnings(as.numeric(get_col(df, c("mirnas_n", "mirna_n", "n_mirnas"))))
mirna_names <- get_col(df, c("mirna_names", "mirnas", "mirna_list"))
p_sep       <- suppressWarnings(as.numeric(get_col(df, c("p_value_separate", "p_value", "p_separate"))))
merged_p    <- suppressWarnings(as.numeric(get_col(df, c("merged_p_value", "p_merged", "combined_p_value"))))
merged_fdr  <- suppressWarnings(as.numeric(get_col(df, c("merged_fdr", "fdr_merged", "adjusted_p_value", "adj_p"))))

# If miRNA count missing or length-0, derive from names
if (is.null(mirnas_n) || length(mirnas_n) == 0) {
  mirnas_n <- derive_miRNA_n(mirna_names, nrow(df))
}


# --- after you've read and cleaned df_raw -> df ---

# Parse "Target Genes (n)" which may be a comma/semicolon-separated list like "3,3,4"
targets_col_raw <- get_col(df, c("target_genes_n", "targets_n", "target_genes"))

parse_counts <- function(x) {
  if (is.null(x)) return(rep(NA_real_, nrow(df)))
  sapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    vals <- suppressWarnings(as.numeric(unlist(stringr::str_split(s, "[;,]"))))
    # choose an aggregator: max is conservative for "any miRNA hits this many genes"
    if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
  })
}

targets_n_agg <- parse_counts(targets_col_raw)

# Keep your other columns as before...
term_full   <- get_col(df, c("go_biological_process_term_name", "term_name", "go_term"))
term_genes  <- suppressWarnings(as.numeric(get_col(df, c("term_genes", "genes_in_term", "term_genes_n"))))
mirnas_n    <- suppressWarnings(as.numeric(get_col(df, c("mirnas_n", "mirna_n", "n_mirnas"))))
mirna_names <- get_col(df, c("mirna_names", "mirnas", "mirna_list"))
merged_p    <- suppressWarnings(as.numeric(get_col(df, c("merged_p_value", "p_merged", "combined_p_value"))))
merged_fdr  <- suppressWarnings(as.numeric(get_col(df, c("merged_fdr", "fdr_merged", "adjusted_p_value", "adj_p"))))

# Derive miRNA count if missing
if (is.null(mirnas_n) || length(mirnas_n) == 0) {
  mirnas_n <- derive_miRNA_n(mirna_names, nrow(df))
}

plot_df <- tibble::tibble(
  term_full   = as.character(term_full),
  term        = stringr::str_trim(as.character(term_full)),
  term_genes  = term_genes,
  targets_n   = targets_n_agg,                # use parsed/aggregated counts
  mirnas_n    = mirnas_n,
  merged_p    = merged_p,
  merged_fdr  = merged_fdr
) |>
  dplyr::filter(!is.na(merged_fdr), is.finite(merged_fdr), merged_fdr > 0) |>
  dplyr::mutate(
    neg_log10_fdr = -log10(merged_fdr),
    hit_fraction  = ifelse(term_genes > 0, targets_n / term_genes, NA_real_)
  )

# Keep only plottable rows
plotable_df <- plot_df |>
  dplyr::filter(
    !is.na(term_full),
    is.finite(neg_log10_fdr),
    !is.na(targets_n), is.finite(targets_n)
  )

# Choose top 30 by smallest FDR, then by larger targets_n
top_n <- 30
df_top <- plotable_df %>%
  arrange(merged_fdr, desc(targets_n)) %>%   # 1) sort by FDR then targets
  slice_head(n = top_n) %>%
  mutate(
    term_label = stringr::str_wrap(term_full, width = 50),
    # 2) reorder factor explicitly by merged_fdr
    term_label = forcats::fct_reorder(term_label, merged_fdr, .desc = TRUE)
  )


# Primary plot: x = -log10(FDR)  (recommended)
p_sig <- ggplot2::ggplot(df_top, aes(x = neg_log10_fdr, y = term_label)) +
  ggplot2::geom_point(aes(size = targets_n, color = neg_log10_fdr)) +
  ggplot2::scale_size_area(max_size = 10, limits = c(0, NA), name = "Observed genes") +
  ggplot2::scale_x_continuous(
    limits = c(4, NA),
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  ggplot2::labs(
    x = expression(-log[10]~"(p-value)"),
    y = NULL,
    color = expression(-log[10]~"(p-value)"),
    title = "Overrepresented GO terms at 12 months",
    subtitle = paste0("Top ", top_n, " terms by unadjusted p-value")
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title  = ggplot2::element_text(face = "bold")
  )

print(p_sig)
ggplot2::ggsave("GO_dotplot_TDP_miRNA_12m.svg", plot = p_sig, width = 8, height = 7, device = "svg")

################################################################################
# APA 6m
################################################################################

# 1) Read Excel
file  <- "PANTHER_GO-Slim_APA_mRNA_6m.xlsx"
sheet <- 1
df_raw <- read_excel(file, sheet = sheet)
df <- janitor::clean_names(df_raw)

# 2) Build tidy table (match PANTHER headers)
observed_vec <- df %>%
  dplyr::select(starts_with("apa_regulated_genes_")) %>%  # first 'APA-regulated genes (...)' col
  dplyr::pull(1) %>%
  suppressWarnings(as.numeric(.))

plot_df <- tibble::tibble(
  term_full    = df$panther_go_slim_biological_process,
  go_id        = stringr::str_extract(df$panther_go_slim_biological_process, "GO:\\d{7}"),
  term         = stringr::str_trim(stringr::str_remove(df$panther_go_slim_biological_process, "\\s*\\(GO:\\d{7}\\)")),
  observed     = observed_vec,
  expected     = suppressWarnings(as.numeric(df$apa_regulated_genes_expected)),
  fold_enrich  = suppressWarnings(as.numeric(df$apa_regulated_genes_fold_enrichment)),
  fdr          = suppressWarnings(as.numeric(df$apa_regulated_genes_fdr))
) %>%
  dplyr::filter(!is.na(fdr), is.finite(fdr), fdr > 0,
                !is.na(observed), is.finite(observed),
                !is.na(fold_enrich), is.finite(fold_enrich)) %>%
  dplyr::mutate(neg_log10_fdr = -log10(fdr))

# 3) Pick how many terms to show
top_n <- 30
df_top <- plot_df %>%
  dplyr::filter(fold_enrich > 1) %>%                         # shows only overrepresented terms, not underrepresented
  dplyr::arrange(fdr, dplyr::desc(fold_enrich)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::mutate(
    term_label = stringr::str_wrap(term, width = 50),
    term_label = forcats::fct_reorder(term_label, neg_log10_fdr)  # order by significance
  )

# Subtitle that adapts if fewer than top_n rows remain
n_show <- nrow(df_top)
subtitle_txt <- if (n_show < top_n) {
  paste0("All ", n_show, " terms by p-value")
} else {
  paste0("Top ", top_n, " terms by p-value")
}

# 4) Dot plot
p_APA_6m <- ggplot(df_top, aes(x = fold_enrich, y = term_label)) +
  geom_point(aes(size = observed, color = neg_log10_fdr)) +
  scale_size_area(max_size = 10, name = "Observed genes") +
  scale_x_continuous(
    limits = c(1, NA),  # force axis to start at 1; let ggplot set the upper limit
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 6)  # ~6 evenly spaced ticks
  ) +
  labs(
    x = "Fold enrichment",
    y = NULL,
    color = expression(-log[10]~p-value),
    title = "Overrepresented GO terms at 6 months",
    subtitle = subtitle_txt
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

print(p_APA_6m)

# 5) Save (fix object name)
ggsave("GO_dotplot_TDP_APA_6m.svg", plot = p_APA_6m, width = 8, height = 9, device = "svg")
