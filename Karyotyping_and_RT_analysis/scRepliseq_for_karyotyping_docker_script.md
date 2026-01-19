#!/bin/bash
################################################################################
# Single-Cell Break Regions Detection Pipeline
#
# Description:
#   Pipeline for identifying chromosome break regions in single-G1 cells
#   generated using scRepli-seq experiment.
#   Processes sequencing data through QC, alignment, and identifies copy number
#   variations (CNVs) / Karyotype, and chromosome breaks using E-Divisive segmentation.
#
# Authors: Jothivanan Elumalai
# License: MIT
# Date: 2025
#
# Citation:
#   If you use this pipeline, please cite:
#   - Elumalai and Hiratani, in revision
#
# Requirements:
#   - Docker (for containerized environment)
#   - Reference genome with BWA index
#   - Blacklist regions file
#
# Docker Image:
#   GitHub: https://github.com/jothivanan-elumalai/single_cell_genome_fragility
#   Image: screpliseq:v1.5
#
# Pipeline Overview:
#   Step 0: Setup directory structure
#   Step 1: Adapter trimming (Trim Galore)
#   Step 2: Read alignment (BWA) and post-processing
#   Step 3: Generate fragment and bin files (AneuFinder)
#   Step 4: Karyotype analysis using E-Divisive segmentation
#   Step 5: Chromosome breakpoint identification
#   NOTE: Steps 0 to 3 as same for both RT and karyotype analysis
#
# Input Files:
#   - Raw sequencing reads (FASTQ.GZ format)
#     * Single-end OR R1 from paired-end (pairs mostly overlap)
#     * File extension MUST be .fastq.gz (not .fq.gz)
#   - Reference genome FASTA with BWA index
#   - Genome index file (.fai)
#   - Blacklist regions (BED format)
#
# Output Structure:
#   project_dir/
#   ├── fastq/              # Raw FASTQ files (user-provided)
#   ├── fastqc/             # Quality control reports
#   ├── trim_fastq/         # Adapter-trimmed reads
#   ├── bam/                # Aligned reads with statistics
#   └── Aneu_analysis/
#       ├── fragment/       # Fragment-level data (GRanges objects)
#       ├── bins/           # Binned genomic data (100kb bins)
#       └── CNV_analysis/
#           ├── models/     # E-Divisive segmentation models
#           └── breaks/     # Identified breakpoint regions
#
# Usage:
#   1. Build Docker image (first time only)
#   2. Edit "Configuration" section below
#   3. Place FASTQ.GZ files in project_dir/fastq/
#   4. Run pipeline step by step
#   5. Monitor quality metrics at each step
#
# Important Notes:
#   - Use R1 reads only if paired-end sequencing
#   - File extension must be .fastq.gz (not .fq.gz)
#   - Minimum of 1 M uniqely mapped reads with MAPQ >=10 is required
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

## ---------------------------------------------------------------
## Configuration - EDIT THESE SETTINGS
## ---------------------------------------------------------------

# Project settings
Project_dir="screpliseq_project"           # Your project directory name
out_name="SAMPLE"                          # Sample/experiment name

# Docker settings
docker_version="screpliseq:v1.5"           # Docker image version
in_dir="${HOME}"                           # Base directory containing project
mount="/data/"                             # Mount point inside container

# Reference genome files (UPDATE THESE PATHS)
genome_name="hg19"                                       # Genome build
genome_index="path/to/reference_genome/hg19.fa"          # BWA-indexed FASTA
genome_file="path/to/reference_genome/hg19.fa.fai"       # Genome index (.fai)
blacklist="path/to/blacklist/hg19-blacklist.bed"         # Blacklist regions

# Adapter sequences
index_seq="AGATCGGAAGAGC"                                 # Illumina adapter
SEQXE="TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN"              # SEQXE protocol

# Computational resources
threads_trim=8          # Threads for trimming (maximum: 8)
threads=20              # Threads for alignment and analysis

## ---------------------------------------------------------------
## Docker Command Setup
## ---------------------------------------------------------------

docker_command="docker run --rm -it -v ${in_dir}:${mount}:rw ${docker_version}"

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

# Verify Docker image exists with "docker images" command

## ---------------------------------------------------------------
## Step 0: Setup Directory Structure
## ---------------------------------------------------------------

echo "Step 0: Setting up directories"
echo "-----------------------------------------------------------"
${docker_command} Step0_setup_directories.sh ${mount}/${Project_dir}
echo "✓ Directory structure created"
echo ""
echo "IMPORTANT: Copy your FASTQ.GZ files to: ${in_dir}/${Project_dir}/fastq/"
echo "Press Enter when files are ready, or Ctrl+C to exit..."
read -r

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
echo "IMPORTANT: Verify that ≥50% of reads have MAPQ > 30"
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
## Step 4: Karyotype Analysis in 500-Kb (E-Divisive Segmentation)
## ---------------------------------------------------------------

echo "Step 4: Karyotype Analysis (E-Divisive)"
echo "-----------------------------------------------------------"
rscript=/usr/local/bin/util/Step9a_R_G1_karyotype_edivisive.R
out_dir="${in_dir}/${Project_dir}/Aneu_analysis/CNV_analysis/models"
mkdir -p "${out_dir}"

files=($(ls ${in_dir}/${Project_dir}/Aneu_analysis/bins/*.Rdata))
karyo_count=0
for file in "${files[@]}"; do
    binfile=$(basename "${file}")
    echo "  Analyzing: $binfile"
    ${docker_command} Rscript --vanilla $rscript \
        ${mount}/${Project_dir}/Aneu_analysis/bins/${binfile} \
        ${mount}/${Project_dir}/Aneu_analysis/CNV_analysis/models \
        ${genome_name}
    karyo_count=$((karyo_count + 1))
done
echo "✓ Karyotype analysis complete for ${karyo_count} samples"
echo "  Models saved in: ${Project_dir}/Aneu_analysis/CNV_analysis/models/"
echo ""

## ---------------------------------------------------------------
## Step 5: Chromosome Break Region Identification
## ---------------------------------------------------------------

echo "Step 5: Chromosome Break Region Identification"
echo "-----------------------------------------------------------"
rscript=/usr/local/bin/util/Step9b_R_G1_karyotype_edivisive_breaks.R
out_dir="${in_dir}/${Project_dir}/Aneu_analysis/CNV_analysis/breaks"
mkdir -p "${out_dir}"

files=($(ls ${in_dir}/${Project_dir}/Aneu_analysis/CNV_analysis/models/*.Rdata))
break_count=0
for file in "${files[@]}"; do
    modelfile=$(basename "${file}")
    echo "  Detecting breakpoints: $modelfile"
    ${docker_command} Rscript --vanilla $rscript \
        ${mount}/${Project_dir}/Aneu_analysis/CNV_analysis/models/${modelfile} \
        ${mount}/${Project_dir}/Aneu_analysis/CNV_analysis/breaks \
        ${mount}/${Project_dir}/Aneu_analysis/fragment
    break_count=$((break_count + 1))
done
echo "✓ Breakpoint detection complete for ${break_count} samples"
echo "  Results in: ${Project_dir}/Aneu_analysis/CNV_analysis/breaks/"
echo ""

## ---------------------------------------------------------------
## End
## ---------------------------------------------------------------