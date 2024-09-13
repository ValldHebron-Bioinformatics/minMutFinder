
![Static Badge](https://img.shields.io/badge/Version-Pre--Release-blue)    ![Static Badge](https://img.shields.io/badge/License-GPL_V3-green)

---
# minMutFinder

**Author**: Ignasi Prats-Méndez  
**Supervisor**: Alejandra González-Sánchez  
**Institution**: HUVH & VHIR  
**Group**: Servei de Microbiologia - Unitat de Virus Respiratoris  

---

## Overview

Welcome to **minMutFinder**, a bioinformatics tool designed to help you find minority mutations in population variants with precision and accuracy. It goes beyond conventional tools by considering the possibility of multiple nucleotide mutations within the same codon. Additionally, it provides detailed metrics regarding your sequences.

If you’re working on viral population data or any genomic dataset, **minMutFinder** can help you achieve deeper insights into minority mutations with a user-friendly and accurate approach.

---

## Features

- **Advanced Mutation Detection**: Identifies minority mutations while accounting for multiple nucleotide changes within a single codon.
- **Comprehensive Analysis**: Provides detailed metrics and plots for a thorough understanding of your sequences.
- **Customizable**: Tailored to support various versions and annotated mutations for enhanced flexibility.

---

## Prerequisites

Before using **minMutFinder**, ensure the following programs are installed:

### Required Software

- **python3**
- **nextflow** v23.10.1 or higher (Installation instructions can be found [here](https://www.nextflow.io/))
- **trimmomatic** v0.39 (via bioconda)
- **minimap2** v2.26-r1175 (via bioconda)
- **lofreq** v2.1.5 (via bioconda)
- **bcftools** v1.17 (via bioconda)
- **samtools** v1.18 (htslib v1.17)

### Required Python Packages

- `os`, `pandas`, `sys`, `csv`, `gzip`, `shutil`, `matplotlib`, `seaborn`, `Bio`, `re`, `plotly`, `numpy`

---

## Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/yourusername/minMutFinder.git
cd minMutFinder
```

Ensure all dependencies are installed via **conda** or **pip**:

```bash
# Install bioconda dependencies
conda install -c bioconda trimmomatic minimap2 lofreq bcftools samtools

# Install Python dependencies
pip install pandas biopython plotly numpy matplotlib seaborn pysam
```

---

## How to Run minMutFinder

To execute **minMutFinder**, follow these steps:

1. **Navigate to the minMutFinder folder** (optional):  
   ```bash
   cd '$path_to_minMutFinder_folder'
   ```

2. **Run minMutFinder using Nextflow**:  
   ```bash
   nextflow run '$path_to_minMutFinder_folder'/minMutFinder.nf --ref_seq <reference.fasta> --out_path <output_name> --r1 <forward_reads.fastq.gz> --r2 <reverse_reads.fastq.gz> --annotate <mutations.tsv> --syn_muts yes
   ```

---

### Arguments

- **`--ref_seq`**: Path and filename of the reference genome FASTA file (1)(2)
- **`--out_path`**: Output name for the virus column
- **`--r1`**: Path and filename of the forward FASTQ compressed file
- **`--r2`**: Path and filename of the reverse FASTQ compressed file
- **`--annotate`**: Path and filename of the TSV file containing the annotated mutations (3)
- **`--syn_muts`**: "yes" or "no", depending on whether to include synonymous mutations in the output plot (default is "no") (4)

---

### Notes

1. The reference genome must contain the coding sequences (CDS) of the proteins. If there are multiple proteins, they should be separated in the FASTA file.
2. FASTA headers must use underscores (`_`) between words. For example: `>NC_006273_2_UL96`.
3. The annotated mutation file should be tab-separated and contain a column named `mutation` for annotated mutations.
4. Use "yes" to include synonymous mutations in the output plot, or "no" to exclude them (default).

---

## License

This project, **minMutFinder**, is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this software under the terms of this license. For more details, refer to the [LICENSE](./LICENSE) file.

---

## Get in Touch

If you encounter any issues, have feature requests, or need assistance, feel free to reach out:

- **Open an issue** directly in this repository by clicking [here](https://github.com/yourusername/minMutFinder/issues).
- **Email us** at [ignasi.prats@vhir.org](mailto:ignasi.prats@vhir.org).

We're always happy to help!

---

## References

- [Nextflow](https://www.nextflow.io/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Lofreq](https://csb5.github.io/lofreq/)
- [Bcftools](http://samtools.github.io/bcftools/bcftools.html)
- [Samtools](http://www.htslib.org/)
