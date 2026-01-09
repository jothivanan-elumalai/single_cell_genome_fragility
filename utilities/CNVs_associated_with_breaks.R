#!/usr/bin/env Rscript
################################################################################
# CNV-Break Overlap Analysis
#
# Description:
#   Analyzes copy number variation (CNV) patterns at genomic break sites by
#   comparing untreated vs. treated cell samples. Identifies CNV gains and
#   losses that occur specifically at stress-induced break regions, focusing
#   on partial overlaps to detect treatment-induced chromosomal changes.
#
# Author: Jothivanan Elumalai
# License: MIT
# Version: 1.0.0
# Date: 2025
#
# Requirements:
#   - R >= 4.0.0
#   - Conda environment: cnv (or similar with required packages)
#
# Required R packages:
#   - data.table
#   - GenomicRanges
#   - ggplot2
#   - ggpubr
#   - tidyverse
#   - rstatix
#
# Installation:
#   install.packages(c("data.table", "ggplot2", "tidyverse"))
#   if (!require("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#   BiocManager::install(c("GenomicRanges"))
#   install.packages(c("ggpubr", "rstatix"))
#
# Input Files:
#   1. Untreated CNV files: BED.GZ files containing CNV calls from control cells
#      Format: chr, start, end, somy/copy_number, [additional columns]
#      Location: Directory containing *500kb_CNV.bed.gz files
#
#   2. Treated CNV files: BED.GZ files containing CNV calls from treated cells
#      Format: chr, start, end, somy/copy_number, [additional columns]
#      Location: Directory containing *500kb_CNV.bed.gz files
#
#   3. Break regions BED file: Genomic regions with stress-induced breaks
#      Format: chr, start, end, break_category
#      Categories: late_CFS, mid_CFS, late_nonCFS, mid_nonCFS, early_nonCFS, early_CFS
#
# Output Files:
#   - [CELLNAME]_CNVs_partial_overlap_with_stress_induced_breaks.txt:
#     Main results table with CNV calls at break regions
#
#   - [CELLNAME]_Aphbreakscat_CNVtype.pdf:
#     Stacked bar plot showing proportion of gains vs. losses by break category
#
#   - [CELLNAME]_Aphbreakscat_CNVsize_byCNVtype.pdf:
#     Box plots comparing CNV sizes across break categories, faceted by gain/loss
#
# Analysis Workflow:
#   1. Load CNV data from untreated cells (control baseline)
#   2. Load CNV data from treated cells (stress-induced changes)
#   3. Load genomic break regions with category annotations
#   4. Determine major (most frequent) copy number for each break region in controls
#   5. Identify CNVs in treated cells that:
#      - Partially overlap break regions (exclude complete containments)
#      - Differ from the major copy number state
#   6. Classify CNVs as gains or losses relative to baseline
#   7. Generate statistical comparisons and visualizations
#
# Key Features:
#   - Focuses on partial overlaps to detect break-associated CNV changes
#   - Excludes CNV segments that fully contain break regions
#   - Determines baseline copy number from untreated cell population
#   - Classifies CNVs as gains (higher) or losses (lower) relative to baseline
#   - Statistical testing: Chi-squared test and pairwise proportion tests
#   - Wilcoxon rank-sum tests for CNV size comparisons
#
# Usage:
#   1. Activate conda environment: conda activate cnv
#   2. Edit file paths in "Configuration" section below
#   3. Run script: Rscript cnv_break_overlap_analysis.R
#   or source("cnv_break_overlap_analysis.R")
#
# Output Interpretation:
#   - CNV type plot: Shows whether breaks are associated with gains or losses
#   - CNV size plot: Compares size distribution of CNVs across break categories
#   - Numbers on plots: Median values (top) and sample counts (bottom)
#   - Statistical values: P-values from Wilcoxon tests comparing categories
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(rstatix)
})

## ---------------------------------------------------------------
## Configuration - EDIT THESE PATHS
## ---------------------------------------------------------------

# Input directories
untreated_dir <- "path/to/control_cnv_files"    # Directory with untreated CNV BED.GZ files
treated_dir   <- "path/to/treated_cnv_files"    # Directory with treated CNV BED.GZ files
break_bed     <- "path/to/break_regions.bed"    # BED file with break regions and categories

# Output settings
out_dir   <- "./cnv_break_analysis_output"      # Output directory
cellname  <- "CELLNAME"                          # Cell line name for output files

# Create output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------
## Read Untreated (Control) CNV Files
## ---------------------------------------------------------------

untreated_files <- list.files(untreated_dir, pattern = "500kb_CNV.bed.gz", full.names = TRUE)
message("Found ", length(untreated_files), " untreated CNV files.")

untreated_grl <- list()
for (file in untreated_files) {
  message("Reading untreated: ", basename(file))

  # Read header
  con <- gzfile(file, "rt")
  header <- readLines(con, n = 1)
  close(con)
  header_cols <- strsplit(header, "\t")[[1]]

  # Read table
  dt <- tryCatch({
    fread(cmd = paste("gzip -cd", shQuote(file)), skip = 1, header = FALSE)
  }, error = function(e) fread(file, skip = 1, header = FALSE))

  if (length(header_cols) == ncol(dt)) setnames(dt, header_cols)
  colnames(dt)[1:4] <- c("chrom", "start", "end", "somy")

  # Clean somy column (handle various naming conventions)
  somy_col <- grep("somy|copy|cn|state|segmean", tolower(names(dt)), value = TRUE)[1]
  somy_raw <- as.character(dt[[somy_col]])
  somy_clean <- gsub("-somy", "", somy_raw)
  somy_clean <- gsub("[^0-9.]", "", somy_clean)
  somy_clean[somy_clean == ""] <- NA
  dt$somy <- as.numeric(somy_clean)
  dt <- dt[!is.na(dt$somy), ]

  # Create GRanges object
  gr <- GRanges(
    seqnames = dt$chrom,
    ranges = IRanges(start = dt$start + 1L, end = dt$end),
    somy = dt$somy
  )
  mcols(gr)$cell_id <- tools::file_path_sans_ext(basename(file))
  untreated_grl[[basename(file)]] <- gr
}

## ---------------------------------------------------------------
## Read Treated CNV Files
## ---------------------------------------------------------------

treated_files <- list.files(treated_dir, pattern = "500kb_CNV.bed.gz", full.names = TRUE)
message("Found ", length(treated_files), " treated CNV files.")

treated_grl <- list()
for (file in treated_files) {
  message("Reading treated: ", basename(file))

  con <- gzfile(file, "rt")
  header <- readLines(con, n = 1)
  close(con)
  header_cols <- strsplit(header, "\t")[[1]]

  dt <- tryCatch({
    fread(cmd = paste("gzip -cd", shQuote(file)), skip = 1, header = FALSE)
  }, error = function(e) fread(file, skip = 1, header = FALSE))

  if (length(header_cols) == ncol(dt)) setnames(dt, header_cols)
  colnames(dt)[1:4] <- c("chrom", "start", "end", "somy")

  somy_col <- grep("somy|copy|cn|state|segmean", tolower(names(dt)), value = TRUE)[1]
  somy_raw <- as.character(dt[[somy_col]])
  somy_clean <- gsub("-somy", "", somy_raw)
  somy_clean <- gsub("[^0-9.]", "", somy_clean)
  somy_clean[somy_clean == ""] <- NA
  dt$somy <- as.numeric(somy_clean)
  dt <- dt[!is.na(dt$somy), ]

  gr <- GRanges(
    seqnames = dt$chrom,
    ranges = IRanges(start = dt$start + 1L, end = dt$end),
    somy = dt$somy
  )
  mcols(gr)$cell_id <- tools::file_path_sans_ext(basename(file))
  treated_grl[[basename(file)]] <- gr
}

## ---------------------------------------------------------------
## Read Break Regions
## ---------------------------------------------------------------

breaks <- fread(break_bed, header = FALSE)
colnames(breaks)[1:4] <- c("chrom", "start", "end", "bcat")
breaks_gr <- GRanges(
  seqnames = breaks$chrom,
  ranges = IRanges(start = breaks$start + 1L, end = breaks$end),
  bcat = breaks$bcat
)
message("Loaded ", length(breaks_gr), " break regions.")

## ---------------------------------------------------------------
## Determine Major (Baseline) Copy Number from Untreated Cells
## ---------------------------------------------------------------

message("Determining baseline copy numbers from untreated cells...")
major_somy <- numeric(length(breaks_gr))

for (i in seq_along(breaks_gr)) {
  br <- breaks_gr[i]
  somy_vals <- c()

  # Collect all copy numbers overlapping this break region across all untreated cells
  for (gr in untreated_grl) {
    hits <- findOverlaps(br, gr)
    if (length(hits) > 0) {
      somy_vals <- c(somy_vals, gr$somy[subjectHits(hits)])
    }
  }

  somy_vals <- somy_vals[!is.na(somy_vals)]

  # Determine most frequent copy number (baseline)
  if (length(somy_vals) > 0) {
    major_somy[i] <- as.numeric(names(sort(table(somy_vals), decreasing = TRUE)[1]))
  } else {
    major_somy[i] <- NA
  }
}

## ---------------------------------------------------------------
## Analyze CNVs in Treated Cells (Partial Overlaps Only)
## ---------------------------------------------------------------

message("Analyzing CNVs in treated cells...")
results <- data.table()

for (i in seq_along(breaks_gr)) {
  br <- breaks_gr[i]
  diff_somy_tbl <- data.table()

  for (gr in treated_grl) {
    hits <- findOverlaps(br, gr)

    if (length(hits) > 0) {
      segs <- gr[subjectHits(hits)]
      segs <- segs[!is.na(segs$somy)]

      # EXCLUDE full containments (CNV segments that fully contain the break region)
      # Keep only partial overlaps to detect break-specific changes
      keep <- !((start(segs) <= start(br) & end(segs) >= end(br)))
      segs <- segs[keep]

      # Keep only CNVs that differ from baseline copy number
      segs <- segs[segs$somy != major_somy[i]]

      if (length(segs) > 0) {
        seg_df <- data.table(
          somy = segs$somy,
          size = width(segs) / 1000000  # Convert to Mb
        )
        diff_somy_tbl <- rbind(diff_somy_tbl, seg_df)
      }
    }
  }

  # Summarize CNVs for this break region
  if (nrow(diff_somy_tbl) > 0) {
    summary_tbl <- diff_somy_tbl[, .(cnv_size = sum(size)), by = somy]
    summary_tbl[, `:=`(
      chr = as.character(seqnames(br)),
      start = start(br),
      end = end(br),
      bcat = elementMetadata(br)$bcat,
      major_somy = major_somy[i]
    )]
    # Classify as gain or loss relative to baseline
    summary_tbl[, cnv_type := ifelse(somy > major_somy, "gain", "loss")]
    results <- rbind(results, summary_tbl)
  } else {
    # No CNVs detected at this break region
    results <- rbind(
      results,
      data.table(
        chr = as.character(seqnames(br)),
        start = start(br),
        end = end(br),
        bcat = elementMetadata(br)$bcat,
        major_somy = major_somy[i],
        somy = NA,
        cnv_size = NA,
        cnv_type = NA
      )
    )
  }
}

## ---------------------------------------------------------------
## Save Results
## ---------------------------------------------------------------

setcolorder(results, c("chr", "start", "end", "bcat", "major_somy", "somy", "cnv_type", "cnv_size"))
output_file <- file.path(out_dir, paste0(cellname, "_CNVs_partial_overlap_with_stress_induced_breaks.txt"))
fwrite(results, output_file, sep = "\t")
message("Results saved to: ", output_file)

## ---------------------------------------------------------------
## Statistical Analysis and Visualization
## ---------------------------------------------------------------

message("Generating plots...")

# Prepare data for plotting
df <- na.omit(results)

# Set factor levels for proper ordering
df$bcat <- factor(df$bcat,
                  levels = c("late_CFS", "mid_CFS", "late_nonCFS",
                            "mid_nonCFS", "early_nonCFS", "early_CFS"))
df$cnv_type <- factor(df$cnv_type, levels = c("loss", "gain"))

## ---------------------------------------------------------------
## Plot 1: CNV Type Distribution (Stacked Bar Plot)
## ---------------------------------------------------------------

count <- na.omit(df) %>%
  group_by(bcat, cnv_type) %>%
  count(cnv_type)

count$bcat <- factor(count$bcat,
                     levels = c("late_CFS", "mid_CFS", "late_nonCFS",
                               "mid_nonCFS", "early_nonCFS", "early_CFS"))
count$cnv_type <- factor(count$cnv_type, levels = c("loss", "gain"))

pdf(file.path(out_dir, paste0(cellname, "_Aphbreakscat_CNVtype.pdf")))
ggplot(data = count, aes(x = bcat, y = n, fill = cnv_type)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Break Categories", y = "CNV Type (%)") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_pubr(base_size = 18, base_family = "Helvetica",
            border = FALSE, legend = c("bottom")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## Statistical Test: Chi-squared test for independence
message("\n--- Chi-squared Test: CNV Type vs Break Category ---")
tab <- xtabs(n ~ bcat + cnv_type, data = count)
# Exclude early_CFS if present
test_result <- chisq.test(tab[!rownames(tab) %in% c("early_CFS"), ])
print(test_result)

## Pairwise proportion tests
message("\n--- Pairwise Proportion Tests (Loss Frequency) ---")
pairwise_result <- pairwise.prop.test(
  x = tab[!rownames(tab) %in% c("early_CFS"), ][, "loss"],
  n = rowSums(tab[!rownames(tab) %in% c("early_CFS"), ]),
  p.adjust.method = "BH"
)
print(pairwise_result)

## ---------------------------------------------------------------
## Plot 2: CNV Size Distribution (Box Plots)
## ---------------------------------------------------------------

median_stats <- na.omit(df) %>%
  group_by(cnv_type, bcat) %>%
  summarise(median_y = median(cnv_size),
           label_txt = round(median_y, 0),
           .groups = "drop")

count_stats <- na.omit(df) %>%
  group_by(cnv_type, bcat) %>%
  count(bcat)

# Comparisons for statistical testing
comp <- list(c("late_CFS", "late_nonCFS"),
            c("mid_CFS", "mid_nonCFS"))

stats <- na.omit(df) %>%
  group_by(cnv_type) %>%
  wilcox_test(cnv_size ~ bcat, comparisons = comp, paired = FALSE) %>%
  mutate(pval = format.pval(p, digits = 2, eps = 0.001))

pdf(file.path(out_dir, paste0(cellname, "_Aphbreakscat_CNVsize_byCNVtype.pdf")),
    width = 8, height = 7)
ggplot(data = df, aes(x = bcat, y = cnv_size, color = cnv_type)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(colour = "grey", alpha = 0.5) +
  stat_pvalue_manual(stats, label = "pval", tip.length = 0,
                    label.size = 6, y.position = c(85, 95)) +
  geom_text(data = median_stats,
           mapping = aes(x = bcat, y = 110, label = label_txt),
           color = "black", size = 6) +
  geom_text(data = count_stats,
           mapping = aes(x = bcat, y = 120, label = n),
           color = "black", size = 6) +
  scale_color_manual(values = c("blue", "red")) +
  facet_wrap(~cnv_type) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  labs(x = "Break Categories", y = "CNV Size (Mb)") +
  theme_pubr(base_size = 18, base_family = "Helvetica",
            border = FALSE, legend = c("none")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

message("\n✓ Analysis complete! Check output directory for results:")
message("  - ", basename(output_file))
message("  - ", paste0(cellname, "_Aphbreakscat_CNVtype.pdf"))
message("  - ", paste0(cellname, "_Aphbreakscat_CNVsize_byCNVtype.pdf"))

################################################################################
## End
################################################################################