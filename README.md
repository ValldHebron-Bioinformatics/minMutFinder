![Static Badge](https://img.shields.io/badge/Version-Pre--Release-blue)    ![Static Badge](https://img.shields.io/badge/License-GPL_V3-green)

---
Title: minMutFinder

Author: Ignasi Prats-Méndez

Supervisor: Alejandra González-Sánchez

Institution: HUVH && VHIR

Group: Servei de Microbiologia - Unitat de Virus Respiratoris

---       

## Welcome to minMutFinder and thanks for using our minority mutations from population variants finder!

You are going to find minMutFinder useful if you are interested on knowing with certainty and accuracy which minority mutations are found in your reads. This tool goes further than most tools, by taking into account the possibility that two different nucleotidic mutations might be situated in the same codon. Moreover, it provides multiple metrics regarding your sequences.

### First of all, you need to have the following programs installed beforehand:
    - python3
    - trimmomatic v0.39 (bioconda)
    - minimap2 v2.26-r1175 (bioconda)
    - lofreq v2.1.5 (bioconda)
    - bcftools v1.17 (bioconda)

### Also the following python packages:
    - sys, pandas, re, Bio, os, turtle, plotly, numpy, gzip, shutil, pysam, csv, matplotlib, seaborn

### Execution:
    - On versions minMutFinder < v0.9.0 or == PYSH 
        - Write on your terminal the following:
            optional: cd '$path_to_minMutFinder_folder'
            python3 '$path_to_minMutFinder_folder'/minMutFinder.py 'arguments'
    - On versions minMutFinder >= v0.9.0
        - Write on your terminal the following:
            optional: cd '$path_to_minMutFinder_folder'
            nextflow run '$path_to_minMutFinder_folder'/minMutFinder.nf 'arguments'

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


#### References 

    trimmomatic
    seqtk
    lofreq
    bcftools
    samtools
    
