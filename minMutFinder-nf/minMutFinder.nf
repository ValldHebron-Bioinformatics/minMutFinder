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

include { dirCreator; refCheck } from './modules/initialization'
include { fastqProcessing; mapping1; variant_calling; read_depth; mapping2; reads_qc_plot } from './modules/from_fq_to_vcf'
include { byProteinAnalysis; protNames } from './modules/by_protein'
include { waitForOutMuts; vizNoAnnot; annotateAndViz } from './modules/viz_and_annotate'

def versionMessage() 
{
	log.info"""
	 
	MINORITY MUTATION FINDER (minMutFinder) - Version: ${workflow.manifest.version} 
	""".stripIndent()
}

def helpMessage() 
{
	log.info"""

![Version](https://img.shields.io/badge/Version-1.0.0-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)

# **minMutFinder**

**Author**: [Ignasi Prats-Méndez](mailto:ignasi.prats@vhir.org)  
**Supervisor**: Alejandra González-Sánchez  
**Institution**: HUVH & VHIR  
**Group**: Servei de Microbiologia - Unitat de Virus Respiratoris  

---

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

**minMutFinder** is a bioinformatics tool designed to help you identify minority mutations in 
population variants with precision and accuracy. Unlike other tools, **minMutFinder** considers 
the possibility of multiple nucleotide mutations within the same codon and provides comprehensive 
metrics for your sequences.

---

## 🔍 Features

- 🧬 **Advanced Mutation Detection**: Identifies minority mutations while accounting for multiple 
nucleotide changes within a single codon.
- 📊 **Comprehensive Analysis**: Provides detailed metrics and plots for a thorough understanding 
of your sequences.
- 🔧 **Customizable**: Supports various versions and annotated mutations for enhanced flexibility.

---

## ⚙️ Arguments

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2)
- `--out_path`: Output name for the virus column
- `--r1`: Path and filename of the forward FASTQ compressed file
- `--r2`: Path and filename of the reverse FASTQ compressed file
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3)
- `--syn_muts`: "yes" or "no", depending on whether to include synonymous mutations in the output 
plot (default is "no") (4)

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

A research paper on **minMutFinder** is currently in progress. In the meantime, please cite this GitHub repository as follows:

> Prats-Méndez I. **minMutFinder**: Minority Mutations Finder. 2024. Available from: https://github.com/ValldHebron-Bioinformatics/minMutFinder

---

## 🔮 Future Work and Limitations

### Current Thresholds:
- **Allele Frequency (AF) ≥ 5%**
- **Read depth per nucleotide position ≥ 20**

### Future Improvements:
- Support for user-provided VCF files along with SAM/BAM files, skipping the initial quality control, mapping, and variant calling steps.
- Currently, **minMutFinder** has been tested with viral datasets. Future developments will expand its capabilities to handle bacterial datasets and larger organisms. There are also plans to explore its potential in oncologic research for cancer datasets.

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
Prints version when asked for
*/

if (params.version || params.v) {
	versionMessage()
	exit 0
}

/**
Prints help when asked for
*/

if (params.help || params.h || params.isEmpty() ) {
	helpMessage()
	exit 0
} else {
    //Creates working dir
    workingpath = params.out_path
    workingdir = file(workingpath)
    if( !workingdir.exists() ) {
    	if( !workingdir.mkdirs() ) 	{
    		exit 1, "Cannot create working directory: $workingpath"
    	} 
    }	
}

// Header log info
log.info """---------------------------------------------
MINORITY MUTATION FINDER (minMutFinder)
---------------------------------------------

Beginning of analysis:

"""

def summary = [:]

summary['Starting time'] = new java.util.Date() 
//Environment
summary['Environment'] = ""
summary['Pipeline Name'] = 'minMutFinder'
summary['Pipeline Version'] = workflow.manifest.version

workflow {
    println "${params}"
    out_paths = dirCreator(params.out_path)
    refCheck(params.ref_seq)
    
    out_fq = fastqProcessing(
                params.r1,
                params.r2,
                params.ref_seq,
                out_paths
    )

    reads_qc_plot(params.r1,params.r2,out_paths, out_fq)
    out_map = mapping1(out_paths, out_fq)

    out_vcf = variant_calling(out_paths, out_map)
    out_map2 = mapping2(out_paths, out_vcf)
    out_depth = read_depth(out_paths, out_map2)
    names = protNames(out_paths, out_vcf)
    out_muts = byProteinAnalysis(out_paths, params.ref_seq, names.flatten(), out_depth)
    out_wait = waitForOutMuts(out_paths, out_vcf, out_muts)
    if ( params.annotate ){
        annotateAndViz(out_paths, params.ref_seq, params.annotate, out_wait, params.syn_muts)
    } else {
        vizNoAnnot(out_paths, params.ref_seq, out_wait, params.syn_muts)
    }

}
