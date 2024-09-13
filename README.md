
![Version](https://img.shields.io/badge/Version-Pre--Release-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)

# minMutFinder

### Author: [Ignasi Prats-Méndez](mailto:ignasi.prats@vhir.org)  
**Supervisor**: Alejandra González-Sánchez  
**Institution**: HUVH && VHIR  
**Group**: Servei de Microbiologia - Unitat de Virus Respiratoris  

---

## Overview

**minMutFinder** is a bioinformatics tool designed to accurately detect minority mutations within population variants from high-throughput sequencing data. It goes beyond most conventional tools by considering the possibility that two different nucleotide mutations may reside in the same codon, providing more precise and reliable results. Additionally, minMutFinder offers a range of useful metrics related to your sequences.

If you're working with viral populations or other genetic datasets and need certainty and accuracy in identifying minority mutations, **minMutFinder** is here to assist you.

---

## Features

- **Advanced Mutation Detection**: Detect minority mutations with accuracy, accounting for multiple mutations in the same codon.
- **Comprehensive Reporting**: Provides detailed metrics and plots for your sequences.
- **Flexibility**: Supports various versions with customizable output and annotated mutation integration.

---

## Requirements

Before running **minMutFinder**, ensure you have the following dependencies installed:

### Software:
- **python3**
- **trimmomatic** v0.39 (via bioconda)
- **minimap2** v2.26-r1175 (via bioconda)
- **lofreq** v2.1.5 (via bioconda)
- **bcftools** v1.17 (via bioconda)

### Python Packages:
- `sys`, `pandas`, `re`, `Bio`, `os`, `turtle`, `plotly`, `numpy`, `gzip`, `shutil`, `pysam`, `csv`, `matplotlib`, `seaborn`

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
conda install -c bioconda trimmomatic minimap2 lofreq bcftools

# Install Python dependencies
pip install pandas biopython plotly numpy matplotlib seaborn pysam
```

---

## Usage

To run **minMutFinder**, follow these guidelines based on the version you are using:

### Arguments

- `-ref-genome`: Path and filename of the reference genome in FASTA format (1)(2)
- `--out-name`: Output name for the virus column
- `--r1`: Path and filename of the forward FASTQ compressed file
- `--r2`: Path and filename of the reverse FASTQ compressed file
- `--annotate`: Path and filename of the TSV file containing annotated mutations (3)
- `--syn_muts`: “yes” or “no” depending on whether to include synonymous mutations in the output plot (default is "no") (4)

### Example Commands

```bash
# For minMutFinder < v0.9.0 or == PYSH
python minMutFinder.py --ref-genome <genome.fasta> --r1 <reads_R1.fastq.gz> --r2 <reads_R2.fastq.gz> --out-name <output_name>

# For minMutFinder >= v0.9.0
python minMutFinder.py --ref-genome <genome.fasta> --r1 <reads_R1.fastq.gz> --r2 <reads_R2.fastq.gz> --out-name <output_name> --annotate <mutations.tsv> --syn_muts <yes/no>
```

### Notes

1. The reference genome must contain individual proteins in separate FASTA entries if applicable.
2. FASTA headers should use underscores (`_`) instead of spaces, e.g., `>NC_006273_2_UL96`.
3. The annotated mutations file should be tab-separated with a column labeled `mutation` for versions >= 0.9.0.
4. Choose "yes" to include synonymous mutations in the final plot.

---

## License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this software under the terms of this license. For more details, refer to the [LICENSE](./LICENSE) file.

---

## Support

We are here to help! If you encounter any issues or have feature requests, please open an [issue](https://github.com/yourusername/minMutFinder/issues) on GitHub or contact us at [ignasi.prats@vhir.org](mailto:ignasi.prats@vhir.org).

---

## References

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Seqtk](https://github.com/lh3/seqtk)
- [Lofreq](https://csb5.github.io/lofreq/)
- [Bcftools](http://samtools.github.io/bcftools/bcftools.html)
- [Samtools](http://www.htslib.org/)

---

## Future Work

We are continuously improving **minMutFinder**. Planned features for upcoming versions include:
- Enhanced visualization tools
- Expanded support for additional mutation types
- Integration with cloud-based platforms for large-scale analysis

Stay tuned for updates!

---

### Final Thoughts

Thank you for using **minMutFinder**! Your contributions and feedback are vital to us. If you have any suggestions for improving this tool, feel free to reach out.
