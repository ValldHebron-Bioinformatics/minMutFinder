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
    - nextflow v23.10.1 or higher
    - trimmomatic v0.39 (bioconda)
    - minimap2 v2.26-r1175 (bioconda)
    - lofreq v2.1.5 (bioconda)
    - bcftools v1.17 (bioconda)
    - samtools v1.18 (htslib v1.17)

### Also the following python packages:
    - os, pandas, sys, csv, gzip, shutil, matplotlib, seaborn, Bio, re, plotly, numpy

### Now you are ready to execute minMutFinder correctly! 


### Execution:
    - Write on your terminal the following:
        optional: cd '$path_to_minMutFinder_folder'
        nextflow run '$path_to_minMutFinder_folder'/minMutFinder.nf 'arguments'

### Arguments:
    - path and filename of the reference genome fasta file (1)(2) = --ref_seq
    - name you want your output to have in the virus column = --out_path
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
<br>
(4)--> "yes" if you want the muations plot to also contain the Synonymous mutations, "no" if you do not want them in the mutations plot

### License

This project, **minMutFinder**, is licensed under the terms of the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

You are free to use, modify, and distribute this software under the terms of the GPL-3.0 license. For more details, please refer to the [LICENSE](./LICENSE) file.


#### References 

    nextflow
    trimmomatic
    lofreq
    bcftools
    samtools
    

