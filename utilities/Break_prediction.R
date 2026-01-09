################################################################################
# Genome-wide Break Prediction Analysis
#
# Description:
#   Logistic regression analysis for predicting genomic break sites using
#   replication timing (RT) landscape features. Analyzes 100kb genomic bins
#   to predict different classes of DNA breaks based on RT derivatives.
#
# Author: Jothivanan Elumalai
# License: MIT
# Version: 1.0.0
#
# Requirements:
#   - R >= 4.0.0
#   - GenomicRanges
#   - dplyr
#   - zoo
#   - nnet
#   - pROC
#   - ggplot2
#   - ggpubr
#   - reshape2
#
# Input Files:
#   1. RT bedGraph file: Replication timing data in bedGraph format
#      Format: chr, start, end, RT_value
#
#   2. Break class BED files: Genomic regions for different break classes
#      Format: chr, start, end, [class_label]
#      Expected classes:
#        - late_CFS (Common Fragile Sites in late S-phase)
#        - mid_CFS (Common Fragile Sites in mid S-phase)
#        - late_nonCFS_delayed (Late S-phase non-CFS, delayed)
#        - late_nonCFS_advanced (Late S-phase non-CFS, advanced)
#        - mid_nonCFS_delayed (Mid S-phase non-CFS, delayed)
#        - mid_nonCFS_advanced (Mid S-phase non-CFS, advanced)
#        - early_nonCFS_delayed (Early S-phase non-CFS, delayed)
#        - early_nonCFS_advanced (Early S-phase non-CFS, advanced)
#
# Output Files:
#   - U2OS_window_features.tsv: Computed features for each genomic bin
#   - U2OS_AUC_model_comparison.tsv: AUC values for all models and classes
#   - U2OS_RT_vs_slope_scatterplot.pdf: RT vs slope visualization
#   - U2OS_RT_vs_curvature_scatterplot.pdf: RT vs curvature visualization
#   - U2OS_AUC_model_comparison.pdf: Model performance comparison
#
# Usage:
#   1. Set working directory and create output folder
#   2. Update file paths in "Read data" section to point to your input files
#   3. Adjust parameters if needed (bin_size, window_bins)
#   4. Run the entire script: source("genome_break_prediction.R")
#
# Analysis Steps:
#   1. Load replication timing and break annotation data
#   2. Overlap genomic bins with break regions
#   3. Calculate RT derivatives (slope and curvature)
#   4. Compute window-based features using rolling windows
#   5. Build logistic regression models with different feature combinations
#   6. Evaluate models using ROC-AUC metrics
#   7. Generate visualizations
#
# Models Evaluated:
#   - RT only: Replication timing alone
#   - slope only: RT slope (rate of change)
#   - curvature only: RT curvature (acceleration)
#   - slope + curvature: Combined derivatives
#   - RT + slope + curvature: Full model
#
################################################################################

## -------------------------
## Setup
## -------------------------
# Load required libraries
library(GenomicRanges)
library(dplyr)
library(zoo)
library(nnet)
library(pROC)
library(ggplot2)
library(ggpubr)
library(reshape2)
options(scipen = 100)

# Create output directory
output_dir <- "./break_prediction_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

## -------------------------
## Parameters
## -------------------------
bin_size <- 100000      # Size of genomic bins (100kb)
window_bins <- 2        # Number of bins on either side for rolling window

## -------------------------
## Read data
## -------------------------
# NOTE: Update these file paths to point to your input data

# Replication timing data
rt <- read.table("path/to/your/RT_data.bedGraph",
                 header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(rt) <- c("chr", "start", "end", "RT")

# Break class annotations
# Update paths to your break annotation files
c1 <- read.table("path/to/late_CFS_breaks.bed",
                 header = FALSE, sep = "\t")
c2 <- read.table("path/to/mid_CFS_breaks.bed",
                 header = FALSE, sep = "\t")
c3d <- read.table("path/to/late_nonCFS_delayed_breaks.bed",
                  header = FALSE, sep = "\t"); c3d$V4 <- "late_nonCFS_delayed"
c3a <- read.table("path/to/late_nonCFS_advanced_breaks.bed",
                  header = FALSE, sep = "\t"); c3a$V4 <- "late_nonCFS_advanced"
c4d <- read.table("path/to/mid_nonCFS_delayed_breaks.bed",
                  header = FALSE, sep = "\t"); c4d$V4 <- "mid_nonCFS_delayed"
c4a <- read.table("path/to/mid_nonCFS_advanced_breaks.bed",
                  header = FALSE, sep = "\t"); c4a$V4 <- "mid_nonCFS_advanced"
c5d <- read.table("path/to/early_nonCFS_delayed_breaks.bed",
                  header = FALSE, sep = "\t"); c5d$V4 <- "early_nonCFS_delayed"
c5a <- read.table("path/to/early_nonCFS_advanced_breaks.bed",
                  header = FALSE, sep = "\t"); c5a$V4 <- "early_nonCFS_advanced"

# Combine all break classes
br <- rbind(c1, c2, c3d, c3a, c4d, c4a, c5d, c5a)
colnames(br) <- c("chr", "start", "end", "break_class")

## -------------------------
## Overlap genomic bins with break regions
## -------------------------
gr_rt <- GRanges(rt$chr, IRanges(rt$start + 1, rt$end))
gr_br <- GRanges(br$chr, IRanges(br$start + 1, br$end),
                 break_class = br$break_class)

hits <- findOverlaps(gr_rt, gr_br)
rt$break_class <- "stable"
rt$break_class[queryHits(hits)] <- br$break_class[subjectHits(hits)]
rt$break_class <- factor(rt$break_class)

## -------------------------
## Calculate RT derivatives
## -------------------------
rt <- rt %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(
    slope = abs((lead(RT) - lag(RT)) / (2 * bin_size)),
    curvature = lead(RT) - 2 * RT + lag(RT)
  ) %>%
  ungroup() %>%
  filter(!is.na(slope), !is.na(curvature))

## -------------------------
## Compute window-based features
## -------------------------
rt_w <- rt %>%
  group_by(chr) %>%
  mutate(
    RT_mean = rollapply(RT, window_bins*2 + 1, mean, fill = NA, align = "center"),
    slope_max = rollapply(abs(slope), window_bins*2 + 1, max, fill = NA, align = "center"),
    curvature_min = rollapply(curvature, window_bins*2 + 1, min, fill = NA, align = "center")
  ) %>%
  ungroup() %>%
  filter(!is.na(RT_mean)) %>%
  mutate(
    RT_z = as.numeric(scale(RT_mean)),
    slope_z = as.numeric(scale(slope_max)),
    curvature_z = as.numeric(scale(curvature_min))
  )

# Save computed features
write.table(
  rt_w,
  file = "computed_window_features.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## -------------------------
## Feature visualization
## -------------------------
# RT vs slope scatter plot
pdf("RT_vs_slope_scatterplot.pdf")
ggplot(rt_w, aes(x = RT_z, y = slope_z)) +
  geom_point(alpha = 0.5) +
  stat_cor(method = "pearson", label.x = min(rt_w$RT_z, na.rm = TRUE),
           label.y = 9, size = 5) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  labs(title = "RT vs Slope Scatter Plot",
       x = "RT (z-score)",
       y = "Slope (z-score)") +
  theme_pubr(base_size = 18, base_family = "Helvetica",
             border = FALSE, legend = c("bottom"))
dev.off()

# RT vs curvature scatter plot
pdf("RT_vs_curvature_scatterplot.pdf")
ggplot(rt_w, aes(x = RT_z, y = curvature_z)) +
  geom_point(alpha = 0.5) +
  stat_cor(method = "pearson", label.x = min(rt_w$RT_z, na.rm = TRUE),
           label.y = -15, size = 5) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(-15, 1), breaks = seq(-15, 0, 2.5)) +
  labs(title = "RT vs Curvature Scatter Plot",
       x = "RT (z-score)",
       y = "Curvature (z-score)") +
  theme_pubr(base_size = 18, base_family = "Helvetica",
             border = FALSE, legend = c("bottom"))
dev.off()

## -------------------------
## Model definitions
## -------------------------
model_list <- list(
  RT = c("RT_z"),
  slope = c("slope_z"),
  curvature = c("curvature_z"),
  slope_curvature = c("slope_z", "curvature_z"),
  RT_slope_curvature = c("RT_z", "slope_z", "curvature_z")
)

## -------------------------
## Build models and evaluate performance
## -------------------------
classes <- levels(rt_w$break_class)

auc_mat <- matrix(NA,
                  nrow = length(classes),
                  ncol = length(model_list),
                  dimnames = list(classes, names(model_list)))

# Loop through models and break classes
for (model_name in names(model_list)) {

  preds <- model_list[[model_name]]

  for (cl in classes) {

    # Create binary outcome variable
    rt_w$tmp <- ifelse(rt_w$break_class == cl, 1, 0)

    # Fit logistic regression model
    m_bin <- glm(
      as.formula(paste("tmp ~", paste(preds, collapse = " + "))),
      family = binomial,
      data = rt_w
    )

    # Predict probabilities and compute ROC-AUC
    prob <- predict(m_bin, type = "response")
    roc_obj <- roc(rt_w$tmp, prob, quiet = TRUE)

    auc_mat[cl, model_name] <- as.numeric(auc(roc_obj))
  }
}

# Convert to data frame and save
auc_df <- as.data.frame(auc_mat)
auc_df$break_class <- rownames(auc_df)

write.table(
  auc_df,
  file = "AUC_model_comparison.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## -------------------------
## Visualize model performance
## -------------------------
# Prepare data for plotting (slope, curvature, slope+curvature models only)
auc_df1 <- auc_df[, c("break_class", "slope", "curvature", "slope_curvature")]
auc_df1 <- reshape2::melt(auc_df1,
                          id.vars = "break_class",
                          variable.name = "model",
                          value.name = "AUC")

# Set factor levels for proper ordering
auc_df1$model <- factor(auc_df1$model,
                        levels = c("slope", "curvature", "slope_curvature"))
auc_df1$break_class <- factor(auc_df1$break_class,
                              levels = c("stable", "late_CFS", "mid_CFS",
                                       "late_nonCFS_advanced", "late_nonCFS_delayed",
                                       "mid_nonCFS_advanced", "mid_nonCFS_delayed",
                                       "early_nonCFS_advanced", "early_nonCFS_delayed"))

# Create AUC comparison plot
pdf("AUC_model_comparison.pdf")
ggplot(auc_df1, aes(x = break_class, y = AUC, color = model)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("blue", "red", "green")) +
  scale_y_continuous(limits = c(0.45, 0.8), breaks = seq(0.45, 0.8, 0.05)) +
  labs(title = "AUC Comparison Across Models",
       x = "Break Class",
       y = "AUC") +
  theme_pubr(base_size = 18, base_family = "Helvetica",
             border = FALSE, legend = c("bottom")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

cat("\nAnalysis complete! Check output directory for results:\n")
cat("  - computed_window_features.tsv\n")
cat("  - AUC_model_comparison.tsv\n")
cat("  - RT_vs_slope_scatterplot.pdf\n")
cat("  - RT_vs_curvature_scatterplot.pdf\n")
cat("  - AUC_model_comparison.pdf\n")

################################################################################
## End
################################################################################