#!/usr/bin/env nextflow

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
    Welcome to minMutFinder and thanks for using our minority mutations from population variants finder!
    You are going to find minMutFinder useful if you are interested on knowing with certainty and accuracy which
    minority mutations are found in your reads.
    This tool goes further than most tools, by taking into account the
    possibility that two different nucleotidic mutations might be situated in the same codon. Moreover,
    it provides multiple metrics regarding your sequences.
        
    First of all, you need to have the following programs installed beforehand:
        - python3
        - trimmomatic v0.39 (bioconda)
        - minimap2 v2.26-r1175 (bioconda)
        - lofreq v2.1.5 (bioconda)

    Also the following python packages:
        - sys, pandas, re, Bio, os, turtle, plotly, numpy, gzip, shutil, pysam, csv, matplotlib, seaborn
    
    
    To do so, do the following:
        ...commands to make minMutFinder.py an executable file
    
    
    Once this is done, execute the following commands in your terminal
        ...commands to add minMutFinder to the computer environment so it can be executed from any directory
    
    
    Now you are ready to execute minMutFinder correctly! 
    
    
    Execution:
        - On versions minMutFinder < v0.9.0
            - Write on your terminal the following:
                optional: cd '$path_to_minMutFinder_folder'
                python3 '$path_to_minMutFinder_folder'/minMutFinder.py 'arguments'
        - On versions minMutFinder >= v0.9.0
            - Write on your terminal the following:
                optional: cd '$path_to_minMutFinder_folder'
                nextflow run '$path_to_minMutFinder_folder'/minMutFinder.nf 'arguments'
    
    Arguments:
    - path and filename of the reference genome fasta file (1)(2) = -ref-genome
    - name you want your output to have in the virus column = --out-name
    - path and filename of the forward fastq compressed file = --r1
    - path and filename of the reverse fastq compressed file = --r2
    - path and filename of the tsv file containing the annotated mutations (3) = --annotate
    - "yes" or "no", depending on the preference (4) = --syn_muts ("no" as default) 

    (1) --> The reference genome must be of the cds of the protein. If there is more than 1 portein in the genome,
            the fasta file must contain separatedly the proteins

    (2)--> The reference genome fasta headers separation between words muts be '_'.
            e.g. '>NC 006273.2 UL96' --> '>NC_006273_2_UL96'

    (3)--> The annotated mutation file must be tab separated, and contain a column named 'mutation' with all the different annotated mutations. 
             This only applies to versions 0.9.0 or higher of minMutFinder.

    (4)--> "yes" if you want the muations plot to also contain the Synonymous mutations, "no" if you do not want them in the mutations plot
    
    
    References 
    
    trimmomatic
    seqtk
    lofreq
    bowtie2
    bcftools
    samtools
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
