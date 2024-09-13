
![Version](https://img.shields.io/badge/Version-Pre--Release-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)  

<p align="center">
  <img src="https://your-logo-url-here.com/logo.png" alt="minMutFinder Logo" width="200">
</p>

---

# minMutFinder

**Author**: [Ignasi Prats-Méndez](mailto:ignasi.prats@vhir.org)  
**Supervisor**: Alejandra González-Sánchez  
**Institution**: HUVH & VHIR  
**Group**: Servei de Microbiologia - Unitat de Virus Respiratoris  

---

## 📜 Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [How to Run](#how-to-run-minMutFinder)
- [Arguments](#arguments)
- [License](#license)
- [Citing minMutFinder](#citing-minmutfinder)
- [Future Work & Limitations](#future-work-and-limitations)
- [Get in Touch](#get-in-touch)
- [References](#references)

---

## 🌟 Overview

**minMutFinder** is a bioinformatics tool designed to help you find minority mutations in population variants with precision and accuracy. It goes beyond conventional tools by considering the possibility of multiple nucleotide mutations within the same codon. Additionally, it provides detailed metrics regarding your sequences.

---

## 🧬 Features

- **Advanced Mutation Detection**: Identifies minority mutations while accounting for multiple nucleotide changes within a single codon.
- **Comprehensive Analysis**: Provides detailed metrics and plots for a thorough understanding of your sequences.
- **Customizable**: Tailored to support various versions and annotated mutations for enhanced flexibility.

---

## ⚙️ Prerequisites

Before using **minMutFinder**, ensure the following programs are installed:

### Required Software
| Software      | Version | Installation |
| ------------- | ------- | ------------ |
| **Nextflow**  | 23.10.1 or higher | [Install](https://www.nextflow.io/) |
| **Python**    | 3.x    |  |
| **Trimmomatic** | 0.39   | [Bioconda](https://bioconda.github.io/) |
| **Minimap2**  | 2.26    | [Bioconda](https://bioconda.github.io/) |
| **Lofreq**    | 2.1.5   | [Bioconda](https://bioconda.github.io/) |
| **Bcftools**  | 1.17    | [Bioconda](https://bioconda.github.io/) |
| **Samtools**  | 1.18    | [Bioconda](https://bioconda.github.io/) |

### Required Python Packages
```
os, pandas, sys, csv, gzip, shutil, matplotlib, seaborn, Bio, re, plotly, numpy
```

---

## 🛠️ Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/yourusername/minMutFinder.git
cd minMutFinder
```

### Step 2: Install Dependencies
```bash
# Install bioconda dependencies
conda install -c bioconda trimmomatic minimap2 lofreq bcftools samtools

# Install Python dependencies
pip install pandas biopython plotly numpy matplotlib seaborn pysam
```

---

## 🚀 How to Run minMutFinder

```bash
nextflow run '$path_to_minMutFinder_folder'/minMutFinder.nf --ref_seq <reference.fasta> --out_path <output_name> --r1 <forward_reads.fastq.gz> --r2 <reverse_reads.fastq.gz> --annotate <mutations.tsv> --syn_muts yes
```

---

## ⚙️ Arguments

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2)
- `--out_path`: Output name for the virus column
- `--r1`: Path and filename of the forward FASTQ compressed file
- `--r2`: Path and filename of the reverse FASTQ compressed file
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3)
- `--syn_muts`: "yes" or "no", depending on whether to include synonymous mutations in the output plot (default is "no") (4)

---

## 📄 License

This project, **minMutFinder**, is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this software under the terms of this license. For more details, refer to the [LICENSE](./LICENSE) file.

---

## 🖊️ Citing minMutFinder

A research paper on **minMutFinder** is currently in progress. In the meantime, please cite this GitHub repository as follows:

> Prats-Méndez I. **minMutFinder**: Minority Mutations Finder. 2024. Available from: https://github.com/yourusername/minMutFinder

Currently, **minMutFinder** has been tested with viral datasets. In future developments, it will be expanded to handle bacterial datasets and larger organisms. Additionally, there are plans to explore its application in oncologic data for cancer research.

---

## 🔮 Future Work and Limitations

Currently, **minMutFinder** works with predefined thresholds of:
- **Allele Frequency (AF) ≥ 5%**
- **Read depth per nucleotide position ≥ 20**

### Future Improvements:
- In future versions, **minMutFinder** will be able to work directly from user-provided VCF files along with SAM/BAM files. This will significantly reduce runtime by skipping the first step of quality control, mapping, and variant calling.

---

## ✉️ Get in Touch

If you encounter any issues, have feature requests, or need assistance, feel free to reach out:

- **Open an issue** directly in this repository by clicking [here](https://github.com/yourusername/minMutFinder/issues).
- **Email us** at [ignasi.prats@vhir.org](mailto:ignasi.prats@vhir.org).

We're always happy to help!

---

## 📚 References

- [Nextflow](https://www.nextflow.io/)
- [Python](https://www.python.org/)
- [Bioconda](https://bioconda.github.io/)
- [Minimap2](https://github.com/lh3/minimap2)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Lofreq](https://csb5.github.io/lofreq/)
- [Bcftools](http://samtools.github.io/bcftools/bcftools.html)
- [Samtools](http://www.htslib.org/)
