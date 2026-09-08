![Version](https://img.shields.io/badge/Version-1.1.1-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)



# **minMutFinder**

## 📜 Table of Contents
- [🎯 Overview](#-overview)
- [🔍 Features](#-features)
- [🛠 Prerequisites](#-prerequisites)
- [📥 Installation](#-installation)
- [🚀 How to Run](#-how-to-run-minmutfinder)
- [⚙️ Arguments](#-arguments)
  - [📝 Notes](#-notes)
- [🧪 Testing & Continuous Integration](#-testing--continuous-integration)
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
| Software        | Version         | Installation |
| --------------- | --------------- | ------------ |
| **Nextflow**    | ≥ 26.04.6       | [Install](https://www.nextflow.io/docs/latest/install.html) |
| **Python**      | 3.10.14         | |
| **libgcc-ng**   | 12 or higher    | [conda-forge](https://conda-forge.org/) |
| **Trimmomatic** | 0.39            | [Bioconda](https://bioconda.github.io/) |
| **Minimap2**    | 2.28            | [Bioconda](https://bioconda.github.io/) |
| **BBmap**       | 39.26           | [Bioconda](https://bioconda.github.io/) |
| **LoFreq**      | 2.1.5           | [Bioconda](https://bioconda.github.io/) |
| **Bcftools**    | 1.22            | [Bioconda](https://bioconda.github.io/) |
| **Samtools**    | 1.22            | [Bioconda](https://bioconda.github.io/) |

### Required Python Packages
| Package        | Version |
| -------------- | ------- |
| **pandas**     | 2.3.3   |
| **numpy**      | 2.2.6   |
| **pysam**      | 0.22.1  |
| **biopython**  | 1.87    |
| **matplotlib** | 3.9.1   |
| **seaborn**    | 0.13.2  |
| **plotly**     | 6.9.0   |

Additional standard library modules used: `os`, `sys`, `csv`, `gzip`, `shutil`, `re`.

---

## 📥 Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/ValldHebron-Bioinformatics/minMutFinder.git
cd minMutFinder
```

### Step 2: Install Dependencies
```bash
# Install conda-forge and bioconda dependencies
conda install -c conda-forge libgcc-ng>=12
conda install -c bioconda \
    minimap2=2.28 \
    bbmap=39.26 \
    samtools=1.22 \
    bcftools=1.22 \
    lofreq=2.1.5 \
    trimmomatic=0.39

# Install Python dependencies
pip install \
    pandas==2.3.3 \
    numpy==2.2.6 \
    pysam==0.22.1 \
    biopython==1.87 \
    matplotlib==3.9.1 \
    seaborn==0.13.2 \
    plotly==6.9.0
```

---

## 🚀 How to Run minMutFinder

```bash
nextflow run minMutFinder.nf --ref_seq <reference.fasta> --out_path <output_name> --r1 <forward_reads.fastq.gz> --r2 <reverse_reads.fastq.gz> --annotate <mutations.tsv> --syn_muts <"yes"/"no">
```

---

## ⚙️ Arguments

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2)(4)
- `--out_path`: Output name for the virus column
- `--r1`: Path and filename of the forward FASTQ compressed file
- `--r2`: Path and filename of the reverse FASTQ compressed file
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3)
- `--syn_muts`: `"yes"` or `"no"`, depending on whether to include synonymous mutations in the output plot (default is `"no"`)
- `--vcf`: Path and filename of the VCF
- `--areads`: Path and filename of the SAM or BAM file
- `--AF`: Number from 0 to 1 with the desired allele frequency threshold (default `0.05`)
- `--depth`: Integer ≥ 0 with the desired read depth threshold per nucleotide position (default `20`)
- `--threads`: Number of threads to be used (default `4`)
- `--SB`: Integer ≥ 0 with the desired strand bias threshold (default `29`)
- `--mapping`: `"minimap2"` or `"bbmap"`, depending on the desired mapping tool. BBmap is highly recommended for reference sequences < 1000 nt (default `"minimap2"`)

### 📝 Notes

  1. The reference genome must contain the coding sequences (CDS) of the proteins. If there are multiple proteins, they should be separated in the FASTA file.
  2. FASTA headers must use underscores (`_`) between words. For example: `>NC_006273_2_UL96`.
  3. The annotated mutation file should be tab-separated and contain a column named `mutation` for annotated mutations.
  4. If FASTA headers share patterns (e.g. `H3N2_PA` and `H3N2_PA_X` share `H3N2_PA`) they must be differentiated (e.g. `H3N2_PA_prot`, `H3N2_PA_X`).

---

## 🧪 Testing & Continuous Integration

**minMutFinder** includes an automated test suite implemented with [nf-test](https://www.nf-test.com/) (v0.9.5) and executed via GitHub Actions on every pull request and release.

### Synthetic Validation Dataset

Pipeline correctness is validated against a controlled synthetic dataset comprising five respiratory virus protocols — IAV H1N1pdm09, IAV H3N2, IAV B Victoria, HRSV-A and HRSV-B — generated with `art_illumina` (MiSeq, 101 bp paired-end, 500× coverage). Each protocol includes three independent replicates. Mutations were inserted at known allele frequencies (1–100%) using a nested population approach, and ground truth is documented in `tests/data/synthetic_mutations.csv`.

### Test Structure

The test suite (`tests/main.nf.test`) comprises 15 tests: 10 setup tests that generate outputs for replicates 2 and 3, and 5 main validation tests (one per virus protocol). Each main test evaluates all three replicates and verifies:

- Correct output file generation (VCF, consensus FASTA, QC metrics, HTML plots)
- Absence of frameshifts, which would indicate a reference mismatch
- Detection of ≥ 5 consensus mutations (AF ≥ 50%) and ≥ 3 minority mutations (5–50%)
- Allele frequencies within the statistically expected range (3σ Binomial interval at 500×)
- Correct annotation of resistance and epitope markers (`Annotated = "yes"`)
- Deterministic reproducibility across the three replicates (0 pp difference)

### Running the Tests Locally

```bash
# Install nf-test
curl -fsSL https://get.nf-test.com | bash

# Run the full test suite
nf-test test tests/main.nf.test
```

### CI/CD Pipeline

The GitHub Actions workflow (`.github/workflows/ci.yml`) provisions a reproducible environment using **Miniforge** and **mamba**, installs all dependencies at the pinned versions listed in the Prerequisites section, and executes the nf-test suite automatically. The workflow is triggered on every pull request and release to ensure that no changes break the validated behaviour of the pipeline.

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
- Make it possible for the user to choose the number of threads, AF, Read depth per nucleotide position.
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
- [LoFreq](https://csb5.github.io/lofreq/)
- [Bcftools](http://samtools.github.io/bcftools/bcftools.html)
- [Samtools](http://www.htslib.org/)
- [nf-test](https://www.nf-test.com/)
