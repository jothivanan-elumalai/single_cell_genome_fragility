#!/bin/bash
set -euo pipefail

# Mapping , post-mapping processing, filtering and generating RPM per 100 kb windows
# Usage: Mapping_to_RPM.sh <fastq_file> <genome> <index> <thread> <mapq> <blacklist> <window_bed_file> <window>

# Check arguments
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <fastq_file> <genome> <index> <thread> <blacklist> <window_bed_file> <window>"
    echo ""
    echo "Arguments:"
    echo "  fastq_file      : Input FASTQ file (gzipped)"
    echo "  genome          : Genome name (e.g., hg19, mm9)"
    echo "  index           : BWA index prefix"
    echo "  thread          : Number of threads"
    echo "  blacklist       : Blacklist regions BED file"
    echo "  window_bed_file : Window BED file for coverage"
    echo "  window          : Window size name (e.g., 100kb, 1Mb)"
    exit 1
fi

# Parse arguments
file=$1
genome=$2
index=$3
THREAD=$4
blacklist=$5
window_bed_file=$6
window=$7

# Extract prefix from filename
prefix=$(basename "$file")
prefix=${prefix%_index_SEQXE_trimmed.fq.gz}

echo "Processing sample: ${prefix}"
echo "Genome: ${genome}"
echo "Threads: ${THREAD}"
echo "MAPQ threshold: ${MAPQ}"

# Define output files
bamfile="bam/${prefix}.${genome}.bam"
cleanbam="bam/${prefix}.${genome}.clean.bam"
sortbam="bam/${prefix}.${genome}.clean_srt.bam"
markdup="bam/${prefix}.${genome}.clean_srt_markdup.bam"
markdup_met="bam/${prefix}.${genome}.clean_srt_markdup.bam.met"
mapq="bam/${prefix}.${genome}.clean_srt_markdup_MAPQ${MAPQ}.bam"
black="bam/${prefix}.${genome}.clean_srt_markdup_MAPQ${MAPQ}_black.bam"

# Step 1: Align reads with BWA
echo "[$(date)] Step 1: Aligning reads with BWA..."
bwa aln -t "${THREAD}" "${index}" "${file}" | \
bwa samse "${index}" - "${file}" | \
samtools view -Sb - > "${bamfile}"

# Step 2: Clean BAM file
echo "[$(date)] Step 2: Cleaning BAM file..."
picard-tools CleanSam I="${bamfile}" O="${cleanbam}"

# Step 3: Sort BAM file
echo "[$(date)] Step 3: Sorting BAM file..."
samtools sort -@ "${THREAD}" -O 'bam' -T "${cleanbam%.*}_srt.bam" "${cleanbam}" > "${sortbam}"
samtools index "${sortbam}"

# Step 4: Mark duplicates
echo "[$(date)] Step 4: Marking duplicates..."
picard-tools MarkDuplicates \
    I="${sortbam}" \
    O="${markdup}" \
    METRICS_FILE="${markdup_met}" \
    REMOVE_DUPLICATES=false

samtools index -@ "${THREAD}" "${markdup}"
samstat "${markdup}"

# Step 5: Filter by MAPQ
echo "[$(date)] Step 5: Filtering by MAPQ ${MAPQ}..."
samtools view -@ "${THREAD}" -b -F 1024 -q 10 "${markdup}" > "${mapq}"

# Step 6: Remove blacklist regions
echo "[$(date)] Step 6: Removing blacklist regions..."
bedtools intersect -v -a "${mapq}" -b "${blacklist}" > "${black}"
samstat "${black}"

# Step 7: Calculate coverage and normalize
echo "[$(date)] Step 7: Calculating coverage and normalizing..."
read_counts=$(samtools view -c "${black}")
rpm=$(echo "scale=5; 1000000 / ${read_counts}" | bc)

echo "Sample_bam_file: ${black}
total_read: ${read_counts}
rpm_factor: ${rpm}" > "${black%.*}_rpminfo.txt"

bedtools bamtobed -i "${black}" | \
bedtools intersect -c -b - -a "${window_bed_file}" > "${black%.*}_${window}_count.bedGraph"

awk -v counts="${read_counts}" '{OFS="\t"}{rpm=1000000/counts;print $1,$2,$3,$4*rpm}' \
    "${black%.*}_${window}_count.bedGraph" > "${black%.*}_${window}_count_rpm.bedGraph"

sort -k1,1 -k2,2n "${black%.*}_${window}_count_rpm.bedGraph" > "${black%.*}_${window}_count_rpm_sort.bedGraph"

# Step 8: Clean up intermediate files
echo "[$(date)] Step 8: Cleaning up intermediate files..."
rm -f "${cleanbam}" "${sortbam}" "${sortbam}.bai" "${mapq}"
rm -f "${black%.*}_${window}_count.bedGraph" "${black%.*}_${window}_count_rpm.bedGraph"
echo "[$(date)] Cleanup complete. Kept: ${bamfile} ${markdup}, ${black}, and RPM bedGraph files"