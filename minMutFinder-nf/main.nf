#!/usr/bin/env nextflow

// This file is part of minMutFinder.
//
// minMutFinder is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// minMutFinder is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
//
// Copyright (C) 2024 Ignasi Prats Méndez

nextflow.enable.dsl = 2

include { dirCreator; refCheck } from './modules/initialization'
include { fastqProcessing;  mapping1_minimap2; mapping1_bbmap; variant_calling_fq; variant_calling_fq_noSB; variant_calling_areads; variant_calling_areads_noSB; mapping2_minimap2; mapping2_bbmap; read_depth; reads_qc_metrics } from './modules/from_fq_to_vcf'
include { inVcf; inVcf_noSB; inAlignedReads } from './modules/input_bam_sam_vcf'
include { byProteinAnalysis; protNames } from './modules/by_protein'
include { waitForOutMuts; vizNoAnnot; annotateAndViz } from './modules/viz_and_annotate'

def versionMessage() {
    log.info """
    MINORITY MUTATION FINDER (minMutFinder) - Version: ${workflow.manifest.version}
    """.stripIndent()
}

def helpMessage() {
    log.info """
    minMutFinder v1.3.0 under GPL-3.0 license

Author: Ignasi Prats-Méndez
Institution: HUVH & VHIR  
Group: Servei de Microbiologia - Unitat de Virus Respiratoris  

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

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2)(4)
- `--out_path`: Output name for the virus column
- `--r1`: Path and filename of the forward FASTQ compressed file
- `--r2`: Path and filename of the reverse FASTQ compressed file
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3)
- `--syn_muts`: "yes" or "no", depending on whether to include synonymous mutations in the output plot (default is "no")

### 📝 Notes

  1. The reference genome must contain the coding sequences (CDS) of the proteins. If there are multiple proteins, they should be separated in the FASTA file.
  2. FASTA headers must use underscores (`_`) between words. For example: `>NC_006273_2_UL96`.
  3. The annotated mutation file should be tab-separated and contain a column named `mutation` for annotated mutations.
  4. If FASTA headers share patterns (e.g. H3N2_PA and H3N2_PA_X share H3N2_PA) they must be differentiated (e.g. H3N2_PA_prot, H3N2_PA_X).

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

"""
}

/**
 * Prints version when asked for
 */
workflow HELP_AND_VERSION_MESSAGES {
    if (params.version || params.v) {
        versionMessage()
        exit 0
    }

    if (params.help || params.h || params.isEmpty()) {
        helpMessage()
        exit 0
    } else {
        // Creates working dir
        workingpath = params.out_path
        workingdir = file(workingpath)
        if (!workingdir.exists()) {
            if (!workingdir.mkdirs()) {
                exit 1, "Cannot create working directory: $workingpath"
            }
        }
    }
}

/**
 * Check if files exist
 */
workflow FILE_CHECK {
    if (params.r1 && !file("${params.r1}").exists()) {
        exit 1, "File ${params.r1} does not exist. Exiting..."
    }

    if (params.r2 && !file("${params.r2}").exists()) {
        exit 1, "File ${params.r2} does not exist. Exiting..."
    }

    if (params.ref_seq && !file("${params.ref_seq}").exists()) {
        exit 1, "File ${params.ref_seq} does not exist. Exiting..."
    }

    if (params.vcf && !file("${params.vcf}").exists()) {
        exit 1, "File ${params.vcf} does not exist. Exiting..."
    }

    if (params.areads && !file("${params.areads}").exists()) {
        exit 1, "File ${params.areads} does not exist. Exiting..."
    }

    if (params.annotate && !file("${params.annotate}").exists()) {
        exit 1, "File ${params.annotate} does not exist. Exiting..."
    }
}

/**
 * Log information about the workflow
 */
workflow LOG_INFO {
    log.info """---------------------------------------------
    MINORITY MUTATIONS FINDER (minMutFinder)
    ---------------------------------------------

    Beginning of analysis:
    """

    def summary = [:]
    summary['Starting time'] = new java.util.Date()
    summary['Environment'] = ""
    summary['Pipeline Name'] = 'minMutFinder'
    summary['Pipeline Version'] = workflow.manifest.version
}

/**
 * Main workflow
 */
workflow MAIN_WORKFLOW {
    HELP_AND_VERSION_MESSAGES()
    FILE_CHECK()
    LOG_INFO()

    dirCreator(params.out_path)
    refs = refCheck(params.out_path, file(params.ref_seq))
    ref_only = refs.map { it[0] }

    if (file("${params.areads}").exists()) {
        out_map = inAlignedReads(file(params.areads), params.out_path)
        if (file("${params.vcf}").exists()) {
            out_vcf = inVcf_noSB(file(params.vcf), params.out_path, out_map, params.AF, params.depth, ref_only)
            // out_vcf = inVcf(file(params.vcf), params.out_path, out_map, params.AF, params.depth, params.SB, ref_only)
        } else {
            // out_vcf = variant_calling_areads(params.out_path, out_map, params.AF, params.depth, params.SB, ref_only, params.threads)
            out_vcf = variant_calling_areads_noSB(params.out_path, out_map, params.AF, params.depth, ref_only, params.threads)
        }
        out_depth = false
    } else {
        out_fq = fastqProcessing(file(params.r1), file(params.r2), params.out_path)
        if (params.mapping == "minimap2") {
            out_map = mapping1_minimap2(params.out_path, out_fq, ref_only)
        } else if (params.mapping == "bbmap") {
            out_map = mapping1_bbmap(params.out_path, out_fq, ref_only)
        } else {
            log.error "Mapping method not recognized. Please use 'minimap2' or 'bbmap'."
            exit "Mapping method not recognized. Please use 'minimap2' or 'bbmap'."
            exit 1
        }
        // out_vcf = variant_calling_fq(params.out_path, out_map, params.AF, params.depth, params.SB, ref_only, params.threads)
        out_vcf = variant_calling_fq_noSB(params.out_path, out_map, params.AF, params.depth, ref_only, params.threads)
        if (params.mapping == "minimap2") {
            out_map2 = mapping2_minimap2(params.out_path, out_vcf, out_fq)
        } else if (params.mapping == "bbmap") {
            out_map2 = mapping2_bbmap(params.out_path, out_vcf, out_fq)
        } else {
            log.error "Mapping method not recognized. Please use 'minimap2' or 'bbmap'."
            exit "Mapping method not recognized. Please use 'minimap2' or 'bbmap'."
            exit 1
        }
        out_depth = read_depth(params.out_path, out_map2)
    }
    samfile_only = out_vcf.map { it[2] }
    names = protNames(params.out_path, out_vcf)
    out_muts = byProteinAnalysis(params.out_path, names.flatten(), params.AF, params.depth, out_depth, samfile_only, ref_only)
    stopper_only = out_muts.map { it[0] }
    out_wait = waitForOutMuts(params.out_path, stopper_only)

    if (params.annotate) {
        annotateAndViz(params.out_path, file(params.annotate), out_wait, params.syn_muts, ref_only)
    } else {
        vizNoAnnot(params.out_path, out_wait, params.syn_muts, ref_only)
    }

    if (file("${params.r1}").exists() && file("${params.r2}").exists()) {
        reads_qc_metrics(file(params.r1), file(params.r2), params.out_path, out_fq, out_wait)
    }

}

workflow {
    MAIN_WORKFLOW()
}
