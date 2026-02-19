# Single-Cell Genome Fragility

Pipelines used in Elumalai & Hiratani, *in revision*

## Overview

This project provides tools and pipelines for:
- **Karyotyping and RT analysis**:
  - Identifying copy number variations (CNVs) and chromosome break regions in single G1 cells using E-Divisive segmentation
  - Measuring DNA RT at using single cells, small-scale population and E/L BrdU-IP method
- **Hi-C analysis**:
  - Scripts used in Elumalai and Hiratani, in revision, for Hi-C analysis using [HiC-Pro](https://github.com/nservant/HiC-Pro)
- **Utilities**: Other scripts used in Elumalai & Hiratani, *in revision*

## Documentation (start with these files)

- [scRepli-seq Karyotyping Pipeline](Karyotyping_and_RT_analysis/scRepliseq_for_karyotyping_docker_script.md) - Complete guide for karyotyping analysis
- [scRepli-seq RT Pipeline](Karyotyping_and_RT_analysis/scRepliseq_for_RT_docker_script.md) - Complete guide for copy-number based RT analysis
- [E/L BrdU-IP Repli-seq RT Pipeline](Karyotyping_and_RT_analysis/EL_BrdUIP_repliseq_RT_docker_script.md) - Complete guide for Early/Late BrdU-IP RT analysis

## Input Requirements for scRepli-seq analysis

### Softwares
- Docker (recommended) or Conda

### Sequencing Data
- Raw sequencing reads in **fastq.gz format**
- Single-end OR R1 from paired-end reads (pairs mostly overlap)

### Reference Files
- Reference genome FASTA with BWA index
- Genome index file (`.fai`)
- [Blacklist regions file](https://github.com/Boyle-Lab/Blacklist)


## Updates in scRepli-seq v1.5 include
- Cutadapt replaced with TrimGalore and this enable multithread usage
- Updated packages: bwa (0.7.19), samtools (1.22), picard-tools (2.18.7) and others
- Updated R packages: R (4.4.3) and AneuFinder (1.34.0)
- Added log2repliscore RT and chromosome break region detection analysis
- Updated docker and conda environment files are provided


## Publications & Citation

If you use this pipeline, please cite:

- **Original scRepli-seq method:**
  - Takahashi et al., Nature Genetics 2019
  - Miura et al., Nature Protocols 2020

- **Original E/L BrdU-IP Repli-seq method:**
  - Hiratani et al., PLOS Biology 2008

- **This work (added log2repliscore RT and scRepli-seq for karyotyping pipelines)**
  - Elumalai & Hiratani, *in revision*

## License

MIT License - See [LICENSE](LICENSE) file for details

## Author

**Jothivanan Elumalai**

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Support

For issues, questions, or suggestions, please open an issue on the GitHub repository.
