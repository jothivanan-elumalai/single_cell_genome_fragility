#!/bin/bash
################################################################################
# Single-Cell Replication Timing (scRepli-seq) Analysis Pipeline
#
# Description:
#   Complete pipeline for analyzing single-cell replication timing from
#   sequencing data. Processes FASTQ files through adapter trimming, alignment,
#   quality control, mappability correction, and replication timing analysis. This is an
#   updated script of [scRepliseq v1.4](https://github.com/kuzobuta/scRepliseq-Pipeline)
#
# Updates in scRepliseq v1.5 include:
#   – Cutadapt replaced with TrimGalore that enable multithread usage
#   – Updated packages: bwa (0.7.19), samtools (1.22), picard-tools (2.18.7) and others
#   – Updated R packages: R (4.4.3) and AneuFinder (1.34.0)
#   – Added log2repliscore RT and chromosome break region detection analysis
#   – Updated docker and conda environment files are provided
#   (NOTE: If you do not use docker, install [samstat](https://github.com/timolassmann/samstat) on your own)
#   (NOTE: Below explanation is only for RT analysis using docker)
#
#
# Authors: Jothivanan Elumalai
# License: MIT
# Version: 1.5
# Date: 2025
#
# Citation:
#   If you use this pipeline, please cite:
#   - Takahashi et al., Nature Genetics 2020 (scRepli-seq)
#   - Miura et al., Nature Protocols 2020 (scRepli-seq protocol)
#   - Elumalai & Hiratani, Nature Communications 2026 (log2repliscore RT method)
#
# Requirements:
#   - Docker (for containerized environment)
#
# Input Files:
#   - Raw sequencing reads (FASTQ.GZ format, single-end or R1 from paired-end)
#   - Reference genome FASTA with BWA index
#   - Genome index file (.fai)
#   - Blacklist regions (BED format)
#
# Output Structure:
#   project_dir/
#   ├── fastq/              # Raw FASTQ files
#   ├── fastqc/             # Quality control reports
#   ├── trim_fastq/         # Adapter-trimmed reads
#   ├── bam/                # Aligned reads
#   └── Aneu_analysis/
#       ├── fragment/       # Fragment-level data
#       ├── bins/           # Binned genomic data
#       ├── MAD_score/      # Quality metrics
#       ├── TagDensity/     # Tag density profiles
#       ├── G1_control/     # Control G1 cells
#       ├── HMM/            # HMM binarization results
#       └── Log2_repliscore/ # RT values
#
# Usage:
#   1. Build Docker image (first time only)
#   2. Edit "Configuration" section below
#   3. Place FASTQ files in project_dir/fastq/
#   4. Run pipeline step by step
#   5. Monitor quality metrics at each step
#
# Pipeline Overview:
#   Step 0: Setup directory structure
#   Step 1: Adapter trimming (Trim Galore)
#   Step 2: Read alignment (BWA) and post-processing
#   Step 3: Generate fragment and bin files (AneuFinder)
#   Step 4: Quality control (MAD scores, tag density)
#   Step 5: G1 control preparation for mappability correction
#   Step 6: 2-HMM binarized RT calculation
#   Step 7: Repliscore calculation (% genome replicated)
#   Step 8: log2repliscore RT value computation
#
# Important Notes:
#   - Use R1 reads only if paired-end (pairs mostly overlap)
#   - File extension must be .fastq.gz (not .fq.gz)
#   - Minimum of 1 M uniqely mapped reads with MAPQ >=10 is required
#   - Select G1 cells with most frequent karyotype for control
#   - log2repliscore method allows RT analysis across all cell cycle phases
#
################################################################################

## ---------------------------------------------------------------
## Configuration - EDIT THESE SETTINGS
## ---------------------------------------------------------------

# Project settings
Project_dir="screpliseq_project"  # Your project directory name
out_name="SAMPLE"                  # Sample/experiment name for outputs

# Docker settings
docker_version="screpliseq:v1.5"
in_dir="${HOME}"                   # Base directory containing project
mount="/data/"                     # Mount point inside Docker container

# Reference genome files (UPDATE THESE PATHS)
genome_name="hg19"
genome_index="path/to/reference_genome/hg19.fa"          # BWA-indexed FASTA
genome_file="path/to/reference_genome/hg19.fa.fai"       # Genome index file
blacklist="path/to/blacklist/hg19-blacklist.bed"         # Blacklist regions

# Adapter sequences
index_seq="AGATCGGAAGAGC"                                 # Illumina adapter
SEQXE="TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN"              # SEQXE protocol sequence

# Computational resources
threads_trim=8          # For trimming (max 8)
threads=20              # Adjust based on your system

# Analysis parameters
binsize="binsize_1e+05"  # 100 Kb bins for HMM
somy="2-somy"            # Use "1-somy" for G1 or early S-phase samples
srt_step=1               # 1 for ascending repliscore order, 0 for no ordering

## ---------------------------------------------------------------
## Build Docker Image (First Time Only)
## ---------------------------------------------------------------

echo "=========================================================="
echo "Single-Cell Replication Timing Analysis Pipeline"
echo "=========================================================="
echo ""
echo "Docker Image Setup"
echo "-----------------------------------------------------------"

if [[ ! -d "single_cell_genome_fragility" ]]; then
    echo "Cloning repository and building Docker image..."
    git clone https://github.com/jothivanan-elumalai/single_cell_genome_fragility.git
    cd single_cell_genome_fragility/Karyotyping_and_RT_analysis
    docker build -t ${docker_version} .
    cd ../..
    echo "✓ Docker image built successfully"
else
    echo "✓ Repository already exists, skipping clone"
fi

# Verify screpliseq Docker image exists by typing "docker images" in shell

## ---------------------------------------------------------------
## Docker Command Setup
## ---------------------------------------------------------------

docker_command="docker run --rm -it -v ${in_dir}:${mount}:rw ${docker_version}"

## ---------------------------------------------------------------
## Step 0: Setup Directory Structure
## ---------------------------------------------------------------

echo "Step 0: Setting up directories"
echo "-----------------------------------------------------------"
${docker_command} Step0_setup_directories.sh ${mount}/${Project_dir}
echo "✓ Directory structure created"
echo ""
echo "IMPORTANT: Copy your FASTQ.GZ files to: ${in_dir}/${Project_dir}/fastq/"

## ---------------------------------------------------------------
## Quality Check: FastQC on Raw Reads
## ---------------------------------------------------------------

echo "Running FastQC on raw reads..."
echo "-----------------------------------------------------------"
files=($(ls ${in_dir}/${Project_dir}/fastq/*.fastq.gz))
for file in "${files[@]}"; do
    in_fq=$(basename "$file")
    echo "  Processing: $in_fq"
    ${docker_command} fastqc -o ${mount}/${Project_dir}/fastqc ${mount}/${Project_dir}/fastq/${in_fq}
done
echo "✓ FastQC complete. Check ${Project_dir}/fastqc/ for reports."
echo ""

## ---------------------------------------------------------------
## Step 1: Adapter Trimming
## ---------------------------------------------------------------

echo "Step 1: Adapter Trimming"
echo "-----------------------------------------------------------"
echo "Adapter sequence: ${index_seq}"
echo "SEQXE sequence: ${SEQXE}"
echo "Threads: ${threads_trim}"
echo ""

files=($(ls ${in_dir}/${Project_dir}/fastq/*.fastq.gz))
for file in "${files[@]}"; do
    prefix=$(basename "$file")
    echo "  Trimming: $prefix"
    ${docker_command} Step1_Trim_fastq.sh \
        ${index_seq} \
        ${mount}/${Project_dir}/fastq/${prefix} \
        ${SEQXE} \
        ${mount}/${Project_dir}/trim_fastq \
        ${threads_trim}
done
echo "✓ Adapter trimming complete"
echo ""

## ---------------------------------------------------------------
## Step 2: Read Alignment and Post-Processing
## ---------------------------------------------------------------

echo "Step 2: Read Alignment (BWA)"
echo "-----------------------------------------------------------"
echo "Reference: ${genome_index}"
echo "Threads: ${threads}"
echo ""

files=($(ls ${in_dir}/${Project_dir}/trim_fastq/*SEQXE_trimmed.fq.gz))
for file in "${files[@]}"; do
    FASTQ=$(basename "${file}")
    prefix=${FASTQ%_index_SEQXE_trimmed.fq.gz}
    echo "  Aligning: $prefix"
    ${docker_command} Step2_Mapping.sh \
        ${mount}/${Project_dir}/trim_fastq/${FASTQ} \
        ${mount}/${genome_index} \
        ${genome_name} \
        ${threads} \
        ${mount}/${Project_dir}/bam/${prefix}
done
echo "✓ Alignment complete"
echo ""

## ---------------------------------------------------------------
## Quality Check: Mapping Statistics
## ---------------------------------------------------------------

echo "Running mapping quality checks (samstat)..."
echo "-----------------------------------------------------------"
files=($(ls ${in_dir}/${Project_dir}/bam/*markdup.bam))
for file in "${files[@]}"; do
    BAM=$(basename "${file}")
    echo "  Checking: $BAM"
    ${docker_command} samstat -t ${threads} ${mount}/${Project_dir}/bam/${BAM}
done
echo "✓ Quality check complete"
echo "IMPORTANT: Verify that ≥50% of reads have MAPQ >= 30"
echo ""

## ---------------------------------------------------------------
## Step 3: Generate Fragment and Bin Files
## ---------------------------------------------------------------

echo "Step 3: Generating Fragment and Bin Files (AneuFinder)"
echo "-----------------------------------------------------------"
rscript=/usr/local/bin/util/Step3_R-Aneu-Fragment-bins.R

files=($(ls ${in_dir}/${Project_dir}/bam/*markdup.bam))
for file in "${files[@]}"; do
    bamfile=$(basename "${file}")
    name=${bamfile%.${genome_name}.clean_srt_markdup.bam}
    echo "  Processing: $name"
    ${docker_command} Rscript --vanilla $rscript \
        ${mount}/${Project_dir}/bam/${bamfile} \
        ${mount}/${Project_dir}/Aneu_analysis \
        ${name} \
        ${mount}/${blacklist} \
        ${mount}/${genome_file}
done
echo "✓ Fragment and bin files generated"
echo ""

## ---------------------------------------------------------------
## Step 4a: Quality Control - MAD Scores
## ---------------------------------------------------------------

echo "Step 4a: Quality Control - MAD Score Calculation"
echo "-----------------------------------------------------------"
rscript=/usr/local/bin/util/Step4a_R_MAD_score.R
${docker_command} Rscript --vanilla $rscript \
    ${mount}/${Project_dir}/Aneu_analysis/bins \
    ${mount}/${Project_dir}/Aneu_analysis/MAD_score \
    $out_name
echo "✓ MAD scores calculated"
echo "Check ${Project_dir}/Aneu_analysis/MAD_score/ for quality metrics"
echo ""

## ---------------------------------------------------------------
## Step 4b: Tag Density Profiles (Optional)
## ---------------------------------------------------------------

echo "Step 4b: Tag Density Profiles (Optional)"
echo "-----------------------------------------------------------"
rscript=/usr/local/bin/util/Step4b_TagDensity_w200ks40k.R
${docker_command} Rscript --vanilla ${rscript} \
    ${mount}/${Project_dir}/Aneu_analysis/fragment \
    ${mount}/${Project_dir}/Aneu_analysis/TagDensity \
    ${mount}/${blacklist} \
    ${mount}/${genome_file}
echo "✓ Tag density profiles generated"
echo ""

## ---------------------------------------------------------------
## Step 5a: G1 Control Preparation
## ---------------------------------------------------------------

echo "Step 5a: G1 Control Cell Preparation"
echo "-----------------------------------------------------------"
echo "IMPORTANT: You need to select G1 control cells manually"
echo "1. Run karyotype analysis on candidate G1 cells"
echo "2. Select cells with most frequent karyotype"
echo "3. Update the G1 bin file paths below"
echo ""
echo "Example karyotype check (modify bin_file path):"
echo ""
# Example for two G1 cell - USER MUST UPDATE the bin file for each G1 cell
rscript=/usr/local/bin/util/Step5a_R_G1_karyotype.R

${docker_command} Rscript --vanilla $rscript \
    ${mount}/${Project_dir}/Aneu_analysis/bins/EXAMPLE_G1_cell1_100k_mapq10_blacklist_bin.Rdata \
    ${mount}/${Project_dir}/Aneu_analysis/G1_control

${docker_command} Rscript --vanilla $rscript \
    ${mount}/${Project_dir}/Aneu_analysis/bins/EXAMPLE_G1_cell2_100k_mapq10_blacklist_bin.Rdata \
    ${mount}/${Project_dir}/Aneu_analysis/G1_control

## ---------------------------------------------------------------
## Step 5b: Merge G1 Control Cells
## ---------------------------------------------------------------

echo "Step 5b: Merge Selected G1 Control Cells"
echo "-----------------------------------------------------------"
echo "After karyotype analysis, merge selected G1 fragment files:"
echo ""
# Enter the fragment files for the selected G1 control cells
# Example for two G1 cells - USER MUST UPDATE the fragment files
rscript=/usr/local/bin/util/Step5b_Merge_fragment_Rdata.R
${docker_command} Rscript --vanilla $rscript \
    ${mount}/${Project_dir}/Aneu_analysis/fragment/EXAMPLE_G1_cell1_100k_mapq10_blacklist_fragment.Rdata \
    ${mount}/${Project_dir}/Aneu_analysis/fragment/EXAMPLE_G1_cell2_100k_mapq10_blacklist_fragment.Rdata \
    -o ${mount}/${Project_dir}/Aneu_analysis/G1_control/Merged_control_G1.fragment.Rdata

## ---------------------------------------------------------------
## Step 6: 2-HMM Binarized RT Calculation
## ---------------------------------------------------------------

echo "Step 6: 2-HMM Binarized RT Calculation"
echo "-----------------------------------------------------------"
echo "Computing replicated/unreplicated states using HMM..."
echo "Bin size: ${binsize}"
echo "Somy: ${somy}"
echo ""

rscript=/usr/local/bin/util/Step6_R_Binarization.R

files=($(ls ${in_dir}/${Project_dir}/Aneu_analysis/bins/*.Rdata))
for file in "${files[@]}"; do
    bin_file=$(basename "${file}")
    echo "  Processing: $bin_file"
    ${docker_command} Rscript --vanilla ${rscript} \
        ${mount}/${Project_dir}/Aneu_analysis/bins/${bin_file} \
        ${mount}/${Project_dir}/Aneu_analysis/HMM \
        ${mount}/${Project_dir}/Aneu_analysis/G1_control/Merged_control_G1.fragment.Rdata \
        ${mount}/${genome_file} \
        ${binsize} \
        ${somy} \
        ${mount}/single_cell_genome_fragility/Karyotyping_and_RT_analysis
done
echo "✓ HMM binarization complete"
echo ""

## ---------------------------------------------------------------
## Step 7: Repliscore Calculation
## ---------------------------------------------------------------

echo "Step 7: Repliscore Calculation (% Genome Replicated)"
echo "-----------------------------------------------------------"
rscript=/usr/local/bin/util/Step7_Repliscores.R

${docker_command} Rscript --vanilla $rscript \
    ${mount}/${Project_dir}/Aneu_analysis/HMM/${somy}/Rdata \
    ${mount}/${Project_dir}/Aneu_analysis/HMM/repliscores \
    ${out_name}_${somy} \
    ${srt_step}
echo "✓ Repliscores calculated"
echo ""

## ---------------------------------------------------------------
## Step 8: log2repliscore RT Calculation
## ---------------------------------------------------------------

echo "Step 8: log2repliscore RT Value Computation"
echo "-----------------------------------------------------------"
echo "Computing continuous RT values across cell cycle..."
echo "Window: 200 Kb, Sliding: 40 Kb"
echo ""

rscript=/usr/local/bin/util/Step8_R_log2_repliscore_RT_scores.R

files=($(ls ${in_dir}/${Project_dir}/Aneu_analysis/fragment/*.Rdata))
for file in "${files[@]}"; do
    fragment=$(basename "${file}")
    echo "  Processing: $fragment"
    ${docker_command} Rscript --vanilla $rscript \
        ${mount}/${Project_dir}/Aneu_analysis/fragment/${fragment} \
        ${mount}/${Project_dir}/Aneu_analysis/Log2_repliscore \
        ${mount}/${Project_dir}/Aneu_analysis/G1_control/Merged_control_G1.fragment.Rdata \
        ${mount}/${blacklist} \
        ${mount}/${genome_file} \
        ${mount}/${Project_dir}/Aneu_analysis/HMM/repliscores/${out_name}_${somy}_Repliscores.txt \
        ${mount}/single_cell_genome_fragility/Karyotyping_and_RT_analysis
done
echo "✓ log2repliscore RT values computed"
echo ""

## ---------------------------------------------------------------
## END
## ---------------------------------------------------------------
