# Single-Cell Genome Fragility

A comprehensive bioinformatics pipeline for analyzing karyotype (chromosome copy number), chromosome break regions, and replication timing in single cells using scRepli-seq (single-cell replication timing sequencing) data.

## Overview

This project provides tools and pipelines for:
- **Karyotyping and break detection**: Identifying copy number variations (CNVs) and chromosome break regions in single G1 cells using E-Divisive segmentation
- **Replication Timing (RT) Analysis**: Measuring DNA replication timing at single-cell or small-scale resolution
- **Utilities**: Other scripts used in Elumalai & Hiratani, *in revision*

## Input Requirements for scRepli-seq analysis

### Sequencing Data
- Raw sequencing reads in **fastq.gz format**
- Single-end OR R1 from paired-end reads (pairs mostly overlap)

### Reference Files
- Reference genome FASTA with BWA index
- Genome index file (`.fai`)
- [Blacklist regions file](https://github.com/Boyle-Lab/Blacklist)

## Installation & Usage of scRepli-seq pipeline

### Using Docker (Recommended)

1. Build the Docker image:
```bash
if [[ ! -d "single_cell_genome_fragility" ]]; then
    echo "Cloning repository and building Docker image..."
    git clone https://github.com/jothivanan-elumalai/single_cell_genome_fragility.git
    cd single_cell_genome_fragility/Karyotyping_and_RT_analysis
    docker_version="screpliseq:v1.5"
    docker build -t ${docker_version} .
    cd ../..
    echo "✓ Docker image built successfully"
else
    echo "✓ Repository already exists, skipping clone"
fi
```

### Using Conda

Install the environment using the provided `screpliseq.yml`:
```bash
cd single_cell_genome_fragility/Karyotyping_and_RT_analysis
conda env create -n screpliseq -f screpliseq.yml
conda activate screpliseq
```
If you are using conda, please install [samstat](https://github.com/TimoLassmann/samstat)

## Updates in scRepli-seq v1.5 include
- Cutadapt replaced with TrimGalore and this enable multithread usage
- Updated packages: bwa (0.7.19), samtools (1.22), picard-tools (2.18.7) and others
- Updated R packages: R (4.4.3) and AneuFinder (1.34.0)
- Added log2repliscore RT and chromosome break region detection analysis
- Updated docker and conda environment files are provided

## Documentation

- [scRepli-seq Karyotyping Pipeline](Karyotyping_and_RT_analysis/scRepliseq_for_karyotyping_docker_script.md) - Complete guide for karyotyping analysis
- [scRepli-seq RT Pipeline](Karyotyping_and_RT_analysis/scRepliseq_for_RT_docker_script.md) - Complete guide for copy-number based replication timing analysis

## Publications & Citation

If you use this pipeline, please cite:

- **Original scRepli-seq method:**
  - Takahashi et al., Nature Genetics 2019
  - Miura et al., Nature Protocols 2020

- **This work:**
  - Elumalai & Hiratani, *in revision*

## License

MIT License - See [LICENSE](LICENSE) file for details

## Author

**Jothivanan Elumalai**

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Support

For issues, questions, or suggestions, please open an issue on the GitHub repository.
