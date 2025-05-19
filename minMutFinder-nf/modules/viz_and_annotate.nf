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

process waitForOutMuts {
    errorStrategy 'terminate'
    publishDir "$params.out_path/mutations", mode: 'copy', pattern: "${out_path.baseName}_mutations.csv"
    input:
    path out_path
    file stopper

    output:
    path("${out_path.baseName}_mutations.csv"), optional: true

    script:
    """
    # Check if all mutations have been processed
    if [[ \$(awk 'END {print NR}' ${out_path}/assembly/prot_names.txt) == \$(awk 'END {print NR}' ${stopper}) ]]; then
        rm ${stopper}
        # Merge all mutations files
        for file in "${out_path}/mutations/${out_path.baseName}"_*_mutations.csv; do
            if [[ -f \$file ]]; then
                tail -n +2 \$file >> ${out_path.baseName}_mutations.csv
                rm \$file
            fi
        done
        cat ${out_path.baseName}_mutations.csv

        # Add header to merged mutations file
        printf "SampleID;Protein;Mutation_type;Aa_change;Amino_Acid_Property_Change;Nt_mutation;Mutation_frequency;Mutation_depth\\n\$(cat ${out_path.baseName}_mutations.csv)" > temp && mv temp ${out_path.baseName}_mutations.csv
    fi
    """
}

process vizNoAnnot {
    errorStrategy 'terminate'
    input:
    path out_path
    file muts
    val syn_muts
    file ref_seq

    output:
    tuple path("${out_path}/plots/Mutations_summary.html"), path("${out_path}/plots/qc_metrics_summary.html")

    script:
    """
    # Generate visualization without annotation
    python3 ${params.project_data}/bin/viz_no_annotate.py --out-dir ${out_path} --ref-seq ${ref_seq} --prot-names ${out_path}/assembly/prot_names.txt --syn-muts ${params.syn_muts}
    rm ${out_path}/assembly/prot_names.txt
    """
}

process annotateAndViz {
    errorStrategy 'terminate'
    input:
    path out_path
    file annotate
    file muts
    val syn_muts
    file ref_seq

    output:
    tuple path("${out_path}/plots/Mutations_summary.html"), path("${out_path}/plots/qc_metrics_summary.html")

    script:
    """
    # Generate visualization with annotation
    python3 ${params.project_data}/bin/annotate_and_viz.py --out-dir ${out_path} --annotate ${params.annotate} --ref-seq ${ref_seq} --prot-names ${out_path}/assembly/prot_names.txt --syn-muts ${params.syn_muts}
    rm ${out_path}/assembly/prot_names.txt
    """
}

/*
process clean_output {
    errorStrategy 'terminate'
    input:
    path out_path
    tuple file(muts_plot), file(qc_plot)

    output:
    stdout

    script:
    """
    DIR_SAMPLE=${out_path}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${out_path}/assembly
    VARIANT_CALLING=${out_path}/variant_calling
    # FASTQ=${out_path}/fastq

    # Clean up intermediate files and directories
    rm -r \$ASSEMBLY
    # if [[ -d \$FASTQ ]]; then
    #     rm -r \$FASTQ
    # fi
    rm \$VARIANT_CALLING/samfile.sam \$VARIANT_CALLING/reads_multifasta_*.fasta \$VARIANT_CALLING/samfile_*.tsv \$VARIANT_CALLING/sam_to_loop_*.sam
    rm \$VARIANT_CALLING/"\$SAMPLE"_indelqual_prot* \$VARIANT_CALLING/"\$SAMPLE"_prot_variants*
    if [[ -d \$VARIANT_CALLING ]] && [[ -z "\$(ls -A "\$VARIANT_CALLING")" ]]; then
        rm -rf \$VARIANT_CALLING
    fi
    """
}
*/