#!/bin/bash
################################################################################
# BrdU-IP Replication Timing (Early/Late) Analysis Pipeline
#
# Description:
#   Pipeline for analyzing replication timing using BrdU-IP sequencing of
#   early and late S-phase cells. Processes early/late S-phase samples
#   through QC, alignment, and calculates log2(Early/Late) replication timing
#   ratios across the genome.
#
# Authors: Jothivanan Elumalai
# License: MIT
# Date: 2025
#
# Citation:
#   If you use this pipeline, please cite:
#   - Elumalai and Hiratani, Nature Communications 2026
#   This method was first described in Hiratani et al., PLOS Biology 2008
#
# Requirements:
#   - Docker (for containerized environment)
#   - Reference genome with BWA index
#   - Blacklist regions file
#   - Genomic windows BED file (e.g., 100kb bins)
#
# Docker Image:
#   GitHub: https://github.com/jothivanan-elumalai/single_cell_genome_fragility
#   Image: screpliseq:v1.5
#
# Pipeline Overview:
#   Step 0: Setup directory structure
#   Step 1: Adapter trimming (Trim Galore)
#   Step 2: Read alignment (BWA) and RPM normalization
#   Step 3: log2(Early/Late) RT value calculation
#
# Experimental Design:
#   This pipeline requires paired samples:
#   - Early S-phase: BrdU-IP from early-S replicating cells
#   - Late S-phase: BrdU-IP from late-S replicating cells
#
#   RT is calculated as log2(Early/Late) ratio:
#   - Positive values: Early replication
#   - Negative values: Late replication
#   - Values near 0: Mid S-phase replication
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

## ---------------------------------------------------------------
## Configuration - EDIT THESE SETTINGS
## ---------------------------------------------------------------

# Project settings
Project_dir="EL_BrdUIP_repliseq_project"  # Your project directory name
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
window_bed_file="path/to/reference_genome/hg19.w100kbs100kb.bed"  # Genomic windows

# Window parameters
window="100kb"          # Window size label (must match window_bed_file)

# Adapter sequences
index_seq="AGATCGGAAGAGC"                                 # Illumina adapter
SEQXE="TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN"              # SEQXE protocol sequence

# Computational resources
threads_trim=8          # For trimming (max 8)
threads=20              # Adjust based on your system

## ---------------------------------------------------------------
## Docker Command Setup
## ---------------------------------------------------------------

docker_command="docker run --rm -it -v ${in_dir}:${mount}:rw ${docker_version}"

## ---------------------------------------------------------------
## Build Docker Image (First Time Only)
## ---------------------------------------------------------------

echo "=========================================================="
echo "Early/late BrdU-IP Replication Timing Analysis Pipeline"
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
## Step 2: Read Alignment and RPM Calculation
## ---------------------------------------------------------------

echo "Step 2: Read Alignment (BWA) and RPM Normalization"
echo "-----------------------------------------------------------"
echo "Reference: ${genome_index}"
echo "Threads: ${threads}"
echo "Window size: ${window}"
echo "Threads: ${threads}"
echo ""

files=($(ls ${in_dir}/${Project_dir}/trim_fastq/*SEQXE_trimmed.fq.gz))
for file in "${files[@]}"; do
    FASTQ=$(basename "${file}")
    prefix=${FASTQ%_index_SEQXE_trimmed.fq.gz}
    echo "  Processing: $prefix"
    ${docker_command} Mapping_to_RPM.sh \
        ${mount}/${Project_dir}/trim_fastq/${FASTQ} \
        ${genome_name} \
        ${mount}/${genome_index} \
        ${threads} \
        ${mount}/${blacklist} \
        ${mount}/${window_bed_file} \
        ${window}
done
echo "✓ Alignment and RPM calculation complete"
echo ""

## ---------------------------------------------------------------
## Step 3: log2(Early/late) RT calculation
## ---------------------------------------------------------------
echo "Step 3: log2(Early/late) RT calculation"
echo "-----------------------------------------------------------"
echo "Do this step for each set of early and late-S RPM files"
echo ""

early="path/to/early_S_RPM_file"          # Early-S RPM per 100-kb file
late="path/to/late_S_RPM_file"          # Late-S RPM per 100-kb file
out_dir="${in_dir}/${Project_dir}/EL_BrdU_RT"
mkdir -p "${out_dir}"

${docker_command} bash -c paste ${early} ${late} | \
awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > \
${mount}/${Project_dir}/${out_dir}/Sample_w100kb_log2EL.bedgraph   # change file name

echo "✓ RT calculation done"
echo ""

## ---------------------------------------------------------------
## END
## ---------------------------------------------------------------
