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
library(tidyverse)

################################################################################
# mRNA 6m
################################################################################

# ---- load data ----
df <- read.csv("DESeq2_mRNA_TDP_6m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Make sure numeric columns are parsed as numeric (handles 1.06E-49 style strings)
num_cols <- c("Counts", "LFC", "SE", "Wald", "p-value", "p-adjusted")
df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters (edit as you like) ----
lfc_thresh  <- 0           # absolute log2 fold-change cutoff
padj_thresh <- 0.05        # FDR cutoff

# ---- prepare data for plotting ----
df_clean <- df %>%
  mutate(`p-adjusted` = ifelse(`p-adjusted` == 0, .Machine$double.xmin, `p-adjusted`)) %>%
  filter(!is.na(LFC), !is.na(`p-adjusted`))

plot_df <- df_clean %>%
  mutate(
    neg_log10_padj = -log10(`p-adjusted`),
    sig = case_when(
      (`p-adjusted` < padj_thresh) & (LFC >=  lfc_thresh) ~ "Up",
      (`p-adjusted` < padj_thresh) & (LFC <= -lfc_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- count up/down ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# ---- volcano plot with counts ----
p_mRNA_6m <- ggplot(plot_df, aes(x = LFC, y = neg_log10_padj)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL
  ) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  labs(
    x = expression(log[2]("fold-change")),
    y = expression(-log[10]("FDR")),
    title = "Dysregulated mRNAs at 6 months"
  ) +
  theme_minimal(base_size = 18) +   # set base font size
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  annotate("text",
           x = max(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_padj, na.rm = TRUE),
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +  # adjust label size
  annotate("text",
           x = min(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_padj, na.rm = TRUE),
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold")

print(p_mRNA_6m)
ggsave("volcano_plot_mRNA_TDP_6m.svg", plot = p_mRNA_6m, width = 8, height = 7, device = "svg")

################################################################################
# mRNA 9m
################################################################################

# ---- load data ----
df <- read.csv("DESeq2_mRNA_TDP_9m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Make sure numeric columns are parsed as numeric (handles 1.06E-49 style strings)
num_cols <- c("Counts", "LFC", "SE", "Wald", "p-value", "p-adjusted")
df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters (edit as you like) ----
lfc_thresh  <- 0           # absolute log2 fold-change cutoff
padj_thresh <- 0.05        # FDR cutoff

# ---- prepare data for plotting ----
df_clean <- df %>%
  mutate(`p-adjusted` = ifelse(`p-adjusted` == 0, .Machine$double.xmin, `p-adjusted`)) %>%
  filter(!is.na(LFC), !is.na(`p-adjusted`))

plot_df <- df_clean %>%
  mutate(
    neg_log10_padj = -log10(`p-adjusted`),
    sig = case_when(
      (`p-adjusted` < padj_thresh) & (LFC >=  lfc_thresh) ~ "Up",
      (`p-adjusted` < padj_thresh) & (LFC <= -lfc_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- count up/down ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# ---- volcano plot with counts ----
p_mRNA_9m <- ggplot(plot_df, aes(x = LFC, y = neg_log10_padj)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL
  ) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  labs(
    x = expression(log[2]("fold-change")),
    y = expression(-log[10]("FDR")),
    title = "Dysregulated mRNAs at 9 months"
  ) +
  theme_minimal(base_size = 18) +   # set base font size
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  annotate("text",
           x = max(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_padj, na.rm = TRUE),
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +  # adjust label size
  annotate("text",
           x = min(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_padj, na.rm = TRUE),
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold")

print(p_mRNA_9m)
ggsave("volcano_plot_mRNA_TDP_9m.svg", plot = p_mRNA_9m, width = 8, height = 7, device = "svg")

################################################################################
# mRNA 12m
################################################################################

# ---- load data ----
df <- read.csv("DESeq2_mRNA_TDP_12m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Make sure numeric columns are parsed as numeric (handles 1.06E-49 style strings)
num_cols <- c("Counts", "LFC", "SE", "Wald", "p-value", "p-adjusted")
df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters (edit as you like) ----
lfc_thresh  <- 0           # absolute log2 fold-change cutoff
padj_thresh <- 0.05        # FDR cutoff

# ---- prepare data for plotting ----
df_clean <- df %>%
  mutate(`p-adjusted` = ifelse(`p-adjusted` == 0, .Machine$double.xmin, `p-adjusted`)) %>%
  filter(!is.na(LFC), !is.na(`p-adjusted`))

plot_df <- df_clean %>%
  mutate(
    neg_log10_padj = -log10(`p-adjusted`),
    sig = case_when(
      (`p-adjusted` < padj_thresh) & (LFC >=  lfc_thresh) ~ "Up",
      (`p-adjusted` < padj_thresh) & (LFC <= -lfc_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- count up/down ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# ---- volcano plot with counts ----
p_mRNA_12m <- ggplot(plot_df, aes(x = LFC, y = neg_log10_padj)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL
  ) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  labs(
    x = expression(log[2]("fold-change")),
    y = expression(-log[10]("FDR")),
    title = "Dysregulated mRNAs at 12 months"
  ) +
  theme_minimal(base_size = 18) +   # set base font size
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  annotate("text",
           x = max(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_padj, na.rm = TRUE),
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +  # adjust label size
  annotate("text",
           x = min(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_padj, na.rm = TRUE),
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold")

print(p_mRNA_12m)
ggsave("volcano_plot_mRNA_TDP_12m.svg", plot = p_mRNA_12m, width = 8, height = 7, device = "svg")

################################################################################
# miRNA 6m
################################################################################

# ---- load data ----
df <- read.csv("DESeq2_miRNA_TDP_6m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Make sure numeric columns are parsed as numeric (handles 1.06E-49 style strings)
num_cols <- c("Counts", "LFC", "SE", "Wald", "p-value", "p-adjusted")
df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters (edit as you like) ----
lfc_thresh  <- 0           # absolute log2 fold-change cutoff
pval_thresh <- 0.05        # unadjusted p-value cutoff

# ---- prepare data for plotting ----
df_clean <- df %>%
  mutate(`p-value` = ifelse(`p-value` == 0, .Machine$double.xmin, `p-value`)) %>%
  filter(!is.na(LFC), !is.na(`p-value`))

plot_df <- df_clean %>%
  mutate(
    neg_log10_pval = -log10(`p-value`),
    sig = case_when(
      (`p-value` < pval_thresh) & (LFC >=  lfc_thresh) ~ "Up",
      (`p-value` < pval_thresh) & (LFC <= -lfc_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- count up/down ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# ---- volcano plot with counts ----
p_miRNA_6m <- ggplot(plot_df, aes(x = LFC, y = neg_log10_pval)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL
  ) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
  labs(
    x = expression(log[2]("fold-change")),
    y = expression(-log[10]("p-value")),
    title = "Dysregulated miRNAs at 6 months"
  ) +
  theme_minimal(base_size = 18) +   # set base font size
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  annotate("text",
           x = max(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_pval, na.rm = TRUE),
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +  # adjust label size
  annotate("text",
           x = min(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_pval, na.rm = TRUE),
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold")

print(p_miRNA_6m)
ggsave("volcano_plot_miRNA_TDP_6m.svg", plot = p_miRNA_6m, width = 8, height = 7, device = "svg")

################################################################################
# miRNA 9m
################################################################################

# ---- load data ----
df <- read.csv("DESeq2_miRNA_TDP_9m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Make sure numeric columns are parsed as numeric (handles 1.06E-49 style strings)
num_cols <- c("Counts", "LFC", "SE", "Wald", "p-value", "p-adjusted")
df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters (edit as you like) ----
lfc_thresh  <- 0           # absolute log2 fold-change cutoff
pval_thresh <- 0.05        # unadjusted p-value cutoff

# ---- prepare data for plotting ----
df_clean <- df %>%
  mutate(`p-value` = ifelse(`p-value` == 0, .Machine$double.xmin, `p-value`)) %>%
  filter(!is.na(LFC), !is.na(`p-value`))

plot_df <- df_clean %>%
  mutate(
    neg_log10_pval = -log10(`p-value`),
    sig = case_when(
      (`p-value` < pval_thresh) & (LFC >=  lfc_thresh) ~ "Up",
      (`p-value` < pval_thresh) & (LFC <= -lfc_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- count up/down ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# ---- volcano plot with counts ----
p_miRNA_9m <- ggplot(plot_df, aes(x = LFC, y = neg_log10_pval)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL
  ) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
  labs(
    x = expression(log[2]("fold-change")),
    y = expression(-log[10]("p-value")),
    title = "Dysregulated miRNAs at 9 months"
  ) +
  theme_minimal(base_size = 18) +   # set base font size
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  annotate("text",
           x = max(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_pval, na.rm = TRUE),
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +  # adjust label size
  annotate("text",
           x = min(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_pval, na.rm = TRUE),
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold")

print(p_miRNA_9m)
ggsave("volcano_plot_miRNA_TDP_9m.svg", plot = p_miRNA_9m, width = 8, height = 7, device = "svg")

################################################################################
# miRNA 12m
################################################################################

# ---- load data ----
df <- read.csv("DESeq2_miRNA_TDP_12m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Make sure numeric columns are parsed as numeric (handles 1.06E-49 style strings)
num_cols <- c("Counts", "LFC", "SE", "Wald", "p-value", "p-adjusted")
df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters (edit as you like) ----
lfc_thresh  <- 0           # absolute log2 fold-change cutoff
pval_thresh <- 0.05        # unadjusted p-value cutoff

# ---- prepare data for plotting ----
df_clean <- df %>%
  mutate(`p-value` = ifelse(`p-value` == 0, .Machine$double.xmin, `p-value`)) %>%
  filter(!is.na(LFC), !is.na(`p-value`))

plot_df <- df_clean %>%
  mutate(
    neg_log10_pval = -log10(`p-value`),
    sig = case_when(
      (`p-value` < pval_thresh) & (LFC >=  lfc_thresh) ~ "Up",
      (`p-value` < pval_thresh) & (LFC <= -lfc_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- count up/down ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# ---- volcano plot with counts ----
p_miRNA_12m <- ggplot(plot_df, aes(x = LFC, y = neg_log10_pval)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL
  ) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
  labs(
    x = expression(log[2]("fold-change")),
    y = expression(-log[10]("p-value")),
    title = "Dysregulated miRNAs at 12 months"
  ) +
  theme_minimal(base_size = 18) +   # set base font size
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  annotate("text",
           x = max(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_pval, na.rm = TRUE),
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +  # adjust label size
  annotate("text",
           x = min(plot_df$LFC, na.rm = TRUE),
           y = max(plot_df$neg_log10_pval, na.rm = TRUE),
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold")

print(p_miRNA_12m)
ggsave("volcano_plot_miRNA_TDP_12m.svg", plot = p_miRNA_12m, width = 8, height = 7, device = "svg")

################################################################################
# APA 6m
################################################################################

# ---- load data ----
apa <- read.csv("test_3UTR_TDP_6m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Ensure numeric columns are parsed correctly
num_cols <- c("RED", "pvalue", "p_adj")
present_num_cols <- intersect(num_cols, names(apa))
apa[present_num_cols] <- lapply(apa[present_num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters ----
red_thresh  <- 0       # effect-size cutoff; 0 means any direction if significant
pval_thresh <- 0.05    # significance cutoff using *pvalue* as requested

# ---- prepare data ----
apa_clean <- apa %>%
  mutate(
    pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)  # avoid -log10(0) = Inf
  ) %>%
  filter(!is.na(RED), !is.na(pvalue))

plot_df <- apa_clean %>%
  mutate(
    neg_log10_pval = -log10(pvalue),
    sig = case_when(
      (pvalue < pval_thresh) & (RED >=  red_thresh) ~ "Up",
      (pvalue < pval_thresh) & (RED <= -red_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- counts ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# Safe top y for annotations (with a little padding)
y_top <- max(plot_df$neg_log10_pval[is.finite(plot_df$neg_log10_pval)], na.rm = TRUE)
y_pad <- diff(range(plot_df$neg_log10_pval[is.finite(plot_df$neg_log10_pval)], na.rm = TRUE)) * 0.02
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2

# ---- volcano plot (RED vs -log10(pvalue)) ----
p_APA_6m <- ggplot(plot_df, aes(x = RED, y = neg_log10_pval)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL  # remove legend title
  ) +
  geom_vline(xintercept = c(-red_thresh, red_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
  labs(
    x = "Relative Expression Difference (RED)",
    y = expression(-log[10]("p-value")),
    title = "APA regulation at 6 months"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.margin = margin(12, 12, 12, 12)
  ) +
  annotate("text",
           x = max(plot_df$RED, na.rm = TRUE),
           y = y_top + y_pad,
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +
  annotate("text",
           x = min(plot_df$RED, na.rm = TRUE),
           y = y_top + y_pad,
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold") +
  coord_cartesian(clip = "off")

print(p_APA_6m)
ggsave("volcano_plot_APA_6m.svg", plot = p_APA_6m, width = 8, height = 7, device = "svg")

################################################################################
# APA 9m
################################################################################

# ---- load data ----
apa <- read.csv("test_3UTR_TDP_9m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Ensure numeric columns are parsed correctly
num_cols <- c("RED", "pvalue", "p_adj")
present_num_cols <- intersect(num_cols, names(apa))
apa[present_num_cols] <- lapply(apa[present_num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters ----
red_thresh  <- 0       # effect-size cutoff; 0 means any direction if significant
pval_thresh <- 0.05    # significance cutoff using *pvalue* as requested

# ---- prepare data ----
apa_clean <- apa %>%
  mutate(
    pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)  # avoid -log10(0) = Inf
  ) %>%
  filter(!is.na(RED), !is.na(pvalue))

plot_df <- apa_clean %>%
  mutate(
    neg_log10_pval = -log10(pvalue),
    sig = case_when(
      (pvalue < pval_thresh) & (RED >=  red_thresh) ~ "Up",
      (pvalue < pval_thresh) & (RED <= -red_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- counts ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# Safe top y for annotations (with a little padding)
y_top <- max(plot_df$neg_log10_pval[is.finite(plot_df$neg_log10_pval)], na.rm = TRUE)
y_pad <- diff(range(plot_df$neg_log10_pval[is.finite(plot_df$neg_log10_pval)], na.rm = TRUE)) * 0.02
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2

# ---- volcano plot (RED vs -log10(pvalue)) ----
p_APA_9m <- ggplot(plot_df, aes(x = RED, y = neg_log10_pval)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL  # remove legend title
  ) +
  geom_vline(xintercept = c(-red_thresh, red_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
  labs(
    x = "Relative Expression Difference (RED)",
    y = expression(-log[10]("p-value")),
    title = "APA regulation at 9 months"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.margin = margin(12, 12, 12, 12)
  ) +
  annotate("text",
           x = max(plot_df$RED, na.rm = TRUE),
           y = y_top + y_pad,
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +
  annotate("text",
           x = min(plot_df$RED, na.rm = TRUE),
           y = y_top + y_pad,
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold") +
  coord_cartesian(clip = "off")

print(p_APA_9m)
ggsave("volcano_plot_APA_9m.svg", plot = p_APA_9m, width = 8, height = 7, device = "svg")

################################################################################
# APA 12m
################################################################################

# ---- load data ----
apa <- read.csv("test_3UTR_TDP_12m.csv", sep = ",", header = TRUE, check.names = FALSE)

# Ensure numeric columns are parsed correctly
num_cols <- c("RED", "pvalue", "p_adj")
present_num_cols <- intersect(num_cols, names(apa))
apa[present_num_cols] <- lapply(apa[present_num_cols], function(x) as.numeric(as.character(x)))

# ---- parameters ----
red_thresh  <- 0       # effect-size cutoff; 0 means any direction if significant
pval_thresh <- 0.05    # significance cutoff using *pvalue* as requested

# ---- prepare data ----
apa_clean <- apa %>%
  mutate(
    pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)  # avoid -log10(0) = Inf
  ) %>%
  filter(!is.na(RED), !is.na(pvalue))

plot_df <- apa_clean %>%
  mutate(
    neg_log10_pval = -log10(pvalue),
    sig = case_when(
      (pvalue < pval_thresh) & (RED >=  red_thresh) ~ "Up",
      (pvalue < pval_thresh) & (RED <= -red_thresh) ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---- counts ----
n_up   <- sum(plot_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(plot_df$sig == "Down", na.rm = TRUE)

# Safe top y for annotations (with a little padding)
y_top <- max(plot_df$neg_log10_pval[is.finite(plot_df$neg_log10_pval)], na.rm = TRUE)
y_pad <- diff(range(plot_df$neg_log10_pval[is.finite(plot_df$neg_log10_pval)], na.rm = TRUE)) * 0.02
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2

# ---- volcano plot (RED vs -log10(pvalue)) ----
p_APA_12m <- ggplot(plot_df, aes(x = RED, y = neg_log10_pval)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70"),
    name = NULL  # remove legend title
  ) +
  geom_vline(xintercept = c(-red_thresh, red_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
  labs(
    x = "Relative Expression Difference (RED)",
    y = expression(-log[10]("p-value")),
    title = "APA regulation at 12 months"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.margin = margin(12, 12, 12, 12)
  ) +
  annotate("text",
           x = max(plot_df$RED, na.rm = TRUE),
           y = y_top + y_pad,
           label = paste("Up =", n_up),
           hjust = 1, vjust = 1, color = "#d73027", size = 6, fontface = "bold") +
  annotate("text",
           x = min(plot_df$RED, na.rm = TRUE),
           y = y_top + y_pad,
           label = paste("Down =", n_down),
           hjust = 0, vjust = 1, color = "#4575b4", size = 6, fontface = "bold") +
  coord_cartesian(clip = "off")

print(p_APA_12m)
ggsave("volcano_plot_APA_12m.svg", plot = p_APA_12m, width = 8, height = 7, device = "svg")
