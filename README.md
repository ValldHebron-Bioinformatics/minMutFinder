![Version](https://img.shields.io/badge/Version-Pre--Release-blue)    ![License](https://img.shields.io/badge/License-GPL_V3-green)

# **minMutFinder**

### Author: Ignasi Prats-Méndez  
### Supervisor: Alejandra González-Sánchez  
### Institution: HUVH && VHIR  
### Group: Servei de Microbiologia - Unitat de Virus Respiratoris  

---

## **Introduction**
Welcome to **minMutFinder** – a powerful tool designed to accurately detect and characterise minority mutations within population variants from Fastq files. If you're looking for a tool that ensures precision by accounting for multiple nucleotide changes within the same codon, **minMutFinder** is the solution you need. The tool also provides comprehensive metrics for your sequences.

---

## **Features**
- **Accurate Minority Mutation Detection**: Detect mutations even at low frequencies with high confidence.
- **Codon-Aware Analysis**: Handles multiple nucleotide mutations within the same codon.
- **Comprehensive Metrics**: Offers detailed sequence metrics for thorough analysis.

---

## **Installation Requirements**

### **Prerequisite Programs** (please ensure the following are installed):
- `python3`
- `trimmomatic v0.39` (via Bioconda)
- `minimap2 v2.26-r1175` (via Bioconda)
- `lofreq v2.1.5` (via Bioconda)
- `bcftools v1.17` (via Bioconda)

### **Required Python Packages**:
- `sys`, `pandas`, `re`, `Bio`, `os`, `turtle`, `plotly`, `numpy`, `gzip`, `shutil`, `pysam`, `csv`, `matplotlib`, `seaborn`

---

## **Usage**

### **Execution**

#### For minMutFinder versions `< v0.9.0` or `== PYSH`:
```bash
cd '$path_to_minMutFinder_folder' # optional
python3 '$path_to_minMutFinder_folder'/minMutFinder.py 'arguments'
```
#### For minMutFinder versions `>= v0.9.0`:
```bash
cd '$path_to_minMutFinder_folder' # optional
nextflow run '$path_to_minMutFinder_folder'/minMutFinder.nf 'arguments'
```

### Arguments:
    - path and filename of the reference genome fasta file (1)(2) = -ref-genome
    - name you want your output to have in the virus column = --out-name
    - path and filename of the forward fastq compressed file = --r1
    - path and filename of the reverse fastq compressed file = --r2
    - path and filename of the tsv file containing the annotated mutations (3) = --annotate
    - "yes" or "no", depending on the preference (4) = --syn_muts ("no" as default) 

(1) --> The reference genome must be of the cds of the protein. If there is more than 1 portein in the genome,
        the fasta file must contain separatedly the proteins
<br>
(2)--> The reference genome fasta headers separation between words muts be '\_'.
        e.g. '>NC 006273.2 UL96' --> '>NC\_006273\_2\_UL96'
<br>
(3)--> The annotated mutation file must be tab separated, and contain a column named 'mutation' with all the different annotated mutations. 
         This only applies to versions 0.9.0 or higher or PYSH of minMutFinder.
<br>
(4)--> "yes" if you want the muations plot to also contain the Synonymous mutations, "no" if you do not want them in the mutations plot

### License

This project, **minMutFinder**, is licensed under the terms of the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

You are free to use, modify, and distribute this software under the terms of the GPL-3.0 license. For more details, please refer to the [LICENSE](./LICENSE) file.


### Get in touch

To report a bug, error, or feature request, please open an issue.

For questions, email us at ignasi.prats@vhir.org; we're happy to help!

#### References 

    trimmomatic
    seqtk
    lofreq
    bcftools
    samtools
    
