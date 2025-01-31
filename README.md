
![Static Badge](https://img.shields.io/badge/Version-Pre--Release-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)

# **minMutFinder**

## 📜 Table of Contents
- [🎯 Overview](#-overview)
- [🔍 Features](#-features)
- [🛠 Prerequisites](#-prerequisites)
- [📥 Installation](#-installation)
- [🚀 How to Run](#-how-to-run-minmutfinder)
- [⚙️ Arguments](#-arguments)
  - [📝 Notes](#-notes)
- [🔏 License](#-license)
- [🖊️ Citing minMutFinder](#-citing-minmutfinder)
- [🔮 Future Work & Limitations](#-future-work-and-limitations)
- [✉️ Get in Touch](#-get-in-touch)
- [📚 References](#-references)

---

## 🎯 Overview

**minMutFinder** is a bioinformatics tool designed to help you identify minority mutations in population variants with precision and accuracy. Unlike other tools, **minMutFinder** considers the possibility of multiple nucleotide mutations within the same codon and provides comprehensive metrics for your sequences.

---

## 🔍 Features

- 🧬 **Advanced Mutation Detection**: Identifies minority mutations while accounting for multiple nucleotide changes within a single codon.
- 📊 **Comprehensive Analysis**: Provides detailed metrics and plots for a thorough understanding of your sequences.
- 🔧 **Customizable**: Supports various versions and annotated mutations for enhanced flexibility.

---

## 🛠 Prerequisites

Ensure the following programs are installed:

### Required Software
| Software      | Version | Installation |
| ------------- | ------- | ------------ |
| **Nextflow**  | 23.10.1 or higher | [Install](https://www.nextflow.io/docs/latest/install.html) |
| **Python**    | 3.9    |  |
| **libgcc-ng** | 12 or higher | [conda-forge](https://conda-forge.org/) |
| **Trimmomatic** | 0.39   | [Bioconda](https://bioconda.github.io/) |
| **Minimap2**  | 2.26    | [Bioconda](https://bioconda.github.io/) |
| **Lofreq**    | 2.1.5   | [Bioconda](https://bioconda.github.io/) |
| **Bcftools**  | 1.17 or higher    | [Bioconda](https://bioconda.github.io/) |
| **Samtools**  | 1.17 or higher   | [Bioconda](https://bioconda.github.io/) |

### Required Python Packages:
```bash
os, pandas, sys, csv, gzip, shutil, matplotlib, seaborn, Bio, re, plotly, numpy
```

---

## 📥 Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/ValldHebron-Bioinformatics/minMutFinder.git
cd minMutFinder
```

### Step 2: Install Dependencies
```bash
# Install conda-forge dependencies
conda install -c conda-forge libgcc-ng>=12

# Install bioconda dependencies
conda install -c bioconda minimap2=2.26 htslib=1.17 bcftools=1.17 samtools=1.18 trimmomatic=0.39 lofreq=2.1.5=py39hb7ef6d5_10

# Install Python dependencies
pip install pandas>=2.1.4 Bio>=1.5.9 plotly>=5.18.0 numpy>=1.26.2 pysam>=0.21.0 matplotlib>=3.8.2 seaborn>=0.13.0
```

---

## 🚀 How to Run minMutFinder

```bash
nextflow run minMutFinder.nf --ref_seq <reference.fasta> --out_path <output_name> --r1 <forward_reads.fastq.gz> --r2 <reverse_reads.fastq.gz> --annotate <mutations.tsv> --syn_muts <"yes"/"no">
```

---

## ⚙️ Arguments

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2)
- `--out_path`: Output name for the virus column
- `--r1`: Path and filename of the forward FASTQ compressed file
- `--r2`: Path and filename of the reverse FASTQ compressed file
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3)
- `--syn_muts`: "yes" or "no", depending on whether to include synonymous mutations in the output plot (default is "no") (4)

### 📝 Notes

  1. The reference genome must contain the coding sequences (CDS) of the proteins. If there are multiple proteins, they should be separated in the FASTA file.
  2. FASTA headers must use underscores (`_`) between words. For example: `>NC_006273_2_UL96`.
  3. The annotated mutation file should be tab-separated and contain a column named `mutation` for annotated mutations.
  4. Use "yes" to include synonymous mutations in the output plot, or "no" to exclude them (default).

---

## 🔏 License

This project, **minMutFinder**, is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this software under the terms of this license. For more details, refer to the [LICENSE](./LICENSE) file.

---

## 🖊️ Citing minMutFinder

A research paper on **minMutFinder** is currently in progress. In the meantime, please cite this GitHub repository using the citation provided by GitHub. You can find the official citation by clicking the **"Cite this repository"** button at the top of the repository page or [view the citation file directly](./CITATION.cff).

---

## 🔮 Future Work and Limitations

### Limitations:
- As of now, **minMutFinder** has only been tested on viral sequencing data.
- Only available for Illumina® sequencing data.

### Current Thresholds:
- **Allele Frequency (AF) ≥ 5%**
- **Read depth per nucleotide position ≥ 20**
- **Number of threads = 64**

### Future Improvements:
- At the moment we are working on uploading **minMutFinder** to nf-core, for easier distribution and use.
- Support for user-provided VCF files along with SAM/BAM files, skipping the initial quality control, mapping, and variant calling steps.
- Make it possible for the user to choose the number of threads, AF, Read depth per nucleotide positon.
- Make it possible for the user to choose kind of output desired (complete or just final result files).
- Support for ONT® sequencing data. 

---

## ✉️ Get in Touch

If you encounter any issues, have feature requests, or need assistance, feel free to reach out:

- **Open an issue** directly in this repository by clicking [here](https://github.com/ValldHebron-Bioinformatics/minMutFinder/issues).
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
