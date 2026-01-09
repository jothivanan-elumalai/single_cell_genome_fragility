#!/bin/bash
################################################################################
# Random Control Region Generator
#
# Description:
#   Generates randomized genomic control regions matched to observed break sites
#   using bedtools shuffle. Creates 100 iterations with different random seeds
#   to generate a robust set of control regions, then categorizes them based on
#   replication timing (RT) to match the distribution of real break sites.
#
# Author: Jothivanan Elumalai
# License: MIT
# Version: 1.0.0
# Date: 2025
#
# Requirements:
#   - bedtools >= 2.30.0
#   - awk (standard Unix tool)
#   - bash >= 4.0
#
# Input Files:
#   1. Break regions BED file: Genomic coordinates of observed break sites
#      Format: chr, start, end, [additional columns]
#
#   2. Genome index file (.fai): Chromosome sizes for the reference genome
#      Format: chr, length, [additional columns]
#      Generated with: samtools faidx reference.fa
#
#   3. Blacklist BED file: Regions to exclude from randomization
#      Format: chr, start, end
#      Examples: centromeres, telomeres, gaps, repetitive regions
#
#   4. Replication timing bedGraph: RT values across the genome
#      Format: chr, start, end, RT_value
#      Used for categorizing controls into early/mid/late S-phase
#
# Output Files:
#   - [PREFIX]_random_control_1to100seed.bed: All 100 randomized iterations
#   - [PREFIX]_random_control_EarlyS.bed: Early S-phase controls (RT > 0.25)
#   - [PREFIX]_random_control_MidS.bed: Mid S-phase controls (-0.25 ≤ RT ≤ 0.25)
#   - [PREFIX]_random_control_LateS.bed: Late S-phase controls (RT < -0.25)
#   - [PREFIX]_random_control_all.bed: Combined file with all categories
#
# Analysis Workflow:
#   1. Create 2 Mb windows across the genome
#   2. Exclude blacklisted regions and original break sites
#   3. Shuffle break sites 100 times with different random seeds
#   4. Ensure no overlaps between randomized regions
#   5. Annotate with replication timing values
#   6. Categorize into early/mid/late S-phase based on RT thresholds
#
# RT Classification Thresholds:
#   - Early S-phase: RT > 0.25 (early replicating)
#   - Mid S-phase: -0.25 ≤ RT ≤ 0.25 (intermediate replicating)
#   - Late S-phase: RT < -0.25 (late replicating)
#
# Key Parameters:
#   - Window size: 2 Mb (for defining eligible regions)
#   - Random seeds: 1-100 (for reproducibility and statistical robustness)
#   - Shuffle options: -noOverlapping (prevents overlapping controls)
#                     -chrom (shuffles within same chromosome)
#
# Usage:
#   1. Edit the "Configuration" section below with your file paths
#   2. Adjust RT thresholds if needed
#   3. Make executable: chmod +x generate_random_controls.sh
#   4. Run: ./generate_random_controls.sh
#
# Notes:
#   - Randomization preserves chromosome and region size
#   - Using 100 seeds provides statistical robustness for downstream analyses
#   - Controls are matched to breaks by RT distribution for fair comparisons
#   - The -noOverlapping flag ensures independent control regions
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

## ---------------------------------------------------------------
## Configuration - EDIT THESE PATHS
## ---------------------------------------------------------------

# Output directory
OUTPUT_DIR="./random_control_output"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Input files - UPDATE THESE PATHS
BREAK_REGIONS="path/to/break_regions.bed"           # Observed break sites
GENOME_INDEX="path/to/reference_genome.fa.fai"      # Genome chromosome sizes
BLACKLIST="path/to/blacklist_regions.bed"           # Regions to exclude
RT_BEDGRAPH="path/to/replication_timing.bedGraph"   # Replication timing data

# Output prefix
OUTPUT_PREFIX="random_control"

# Parameters
WINDOW_SIZE=2000000      # 2 Mb windows for defining eligible regions
NUM_SEEDS=100            # Number of randomization iterations
RT_EARLY_THRESHOLD=0.25  # RT threshold for early S-phase
RT_LATE_THRESHOLD=-0.25  # RT threshold for late S-phase

## ---------------------------------------------------------------
## Validate Input Files
## ---------------------------------------------------------------

echo "==================================================="
echo "Random Control Region Generator"
echo "==================================================="
echo ""

for file in "$BREAK_REGIONS" "$GENOME_INDEX" "$BLACKLIST" "$RT_BEDGRAPH"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        echo "Please update file paths in the Configuration section."
        exit 1
    fi
done

echo "✓ All input files found"
echo "  - Break regions: $BREAK_REGIONS"
echo "  - Genome index: $GENOME_INDEX"
echo "  - Blacklist: $BLACKLIST"
echo "  - RT bedGraph: $RT_BEDGRAPH"
echo ""

## ---------------------------------------------------------------
## Generate Random Controls with Multiple Seeds
## ---------------------------------------------------------------

OUTPUT_ALL="${OUTPUT_PREFIX}_1to${NUM_SEEDS}seed.bed"

echo "Generating randomized controls..."
echo "  - Number of iterations: $NUM_SEEDS"
echo "  - Window size: $WINDOW_SIZE bp"
echo "  - Output file: $OUTPUT_ALL"
echo ""

# Remove output file if it exists
[[ -f "$OUTPUT_ALL" ]] && rm "$OUTPUT_ALL"

# Loop through seeds
for seed in $(seq 1 $NUM_SEEDS); do
    if (( seed % 10 == 0 )); then
        echo "  Processing seed $seed/$NUM_SEEDS..."
    fi

    # Generate eligible regions (2Mb windows excluding blacklist and breaks)
    # Then shuffle break sites within these regions
    bedtools makewindows -g "$GENOME_INDEX" -w $WINDOW_SIZE | \
    bedtools intersect -v -a - -b "$BLACKLIST" "$BREAK_REGIONS" | \
    bedtools shuffle \
        -noOverlapping \
        -chrom \
        -incl - \
        -seed $seed \
        -g "$GENOME_INDEX" \
        -i "$BREAK_REGIONS" | \
    awk '{print $1, $2, $3}' OFS='\t' >> "$OUTPUT_ALL"
done

echo "✓ Generated $(wc -l < "$OUTPUT_ALL") randomized regions"
echo ""

## ---------------------------------------------------------------
## Annotate with Replication Timing and Categorize
## ---------------------------------------------------------------

echo "Categorizing controls by replication timing..."

# Early S-phase (RT > 0.25)
OUTPUT_EARLY="${OUTPUT_PREFIX}_EarlyS.bed"
bedtools sort -i "$OUTPUT_ALL" | \
bedtools map -c 4 -o mean -null NA -a - -b "$RT_BEDGRAPH" | \
awk -v thresh=$RT_EARLY_THRESHOLD \
    '{if($4 > thresh) print $1, $2, $3, $5="early_rc"}' OFS='\t' \
    > "$OUTPUT_EARLY"
echo "  - Early S-phase: $(wc -l < "$OUTPUT_EARLY") regions (RT > $RT_EARLY_THRESHOLD)"

# Late S-phase (RT < -0.25)
OUTPUT_LATE="${OUTPUT_PREFIX}_LateS.bed"
bedtools sort -i "$OUTPUT_ALL" | \
bedtools map -c 4 -o mean -null NA -a - -b "$RT_BEDGRAPH" | \
awk -v thresh=$RT_LATE_THRESHOLD \
    '{if($4 < thresh) print $1, $2, $3, $5="late_rc"}' OFS='\t' \
    > "$OUTPUT_LATE"
echo "  - Late S-phase: $(wc -l < "$OUTPUT_LATE") regions (RT < $RT_LATE_THRESHOLD)"

# Mid S-phase (-0.25 <= RT <= 0.25)
OUTPUT_MID="${OUTPUT_PREFIX}_MidS.bed"
bedtools sort -i "$OUTPUT_ALL" | \
bedtools map -c 4 -o mean -null NA -a - -b "$RT_BEDGRAPH" | \
awk -v thresh_low=$RT_LATE_THRESHOLD -v thresh_high=$RT_EARLY_THRESHOLD \
    '{if($4 <= thresh_high && $4 >= thresh_low) print $1, $2, $3, $5="mid_rc"}' OFS='\t' \
    > "$OUTPUT_MID"
echo "  - Mid S-phase: $(wc -l < "$OUTPUT_MID") regions ($RT_LATE_THRESHOLD ≤ RT ≤ $RT_EARLY_THRESHOLD)"

# Combine all categories
OUTPUT_COMBINED="${OUTPUT_PREFIX}_all.bed"
cat "$OUTPUT_EARLY" "$OUTPUT_MID" "$OUTPUT_LATE" > "$OUTPUT_COMBINED"
echo "  - Combined total: $(wc -l < "$OUTPUT_COMBINED") regions"
echo ""

## ---------------------------------------------------------------
## Summary Statistics
## ---------------------------------------------------------------

echo "==================================================="
echo "Summary"
echo "==================================================="
echo ""
echo "Generated Files:"
echo "  1. $OUTPUT_ALL"
echo "     - All randomized controls (100 iterations)"
echo ""
echo "  2. $OUTPUT_EARLY"
echo "     - Early S-phase controls"
echo ""
echo "  3. $OUTPUT_MID"
echo "     - Mid S-phase controls"
echo ""
echo "  4. $OUTPUT_LATE"
echo "     - Late S-phase controls"
echo ""
echo "  5. $OUTPUT_COMBINED"
echo "     - All categorized controls combined"
echo ""
echo "RT Classification Thresholds:"
echo "  - Early: RT > $RT_EARLY_THRESHOLD"
echo "  - Mid: $RT_LATE_THRESHOLD ≤ RT ≤ $RT_EARLY_THRESHOLD"
echo "  - Late: RT < $RT_LATE_THRESHOLD"
echo ""
echo "✓ Analysis complete!"
echo "==================================================="

################################################################################
## End
################################################################################