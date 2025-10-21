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
// Copyright (C) 2024 Respiratory Viruses Unit, Microbiology Department, Vall d’Hebron Hospital Universitari, Vall d’Hebron Institut de Recerca (VHIR), Vall d’Hebron Barcelona Hospital Campus, Passeig Vall d’Hebron 119-129, 08035 Barcelona, Spain

process waitForOutMuts {
    input:
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    file indelqual_prot_bam
    file indelqual_prot_bam_indexed
    file prot_variants_vcf
    file indelqual_prot_sam
    file samfile
    file prot_variants_vcf_gz
    file prot_variants_vcf_gz_indexed
    file vcf_file
    file prot_variants_AF05_vcf_gz
    file prot_variants_AF05_vcf_gz_indexed
    file prot_variants_AF05_fasta
    file mutations
    file qc_metrics
    file depth

    output:
    path("${muts_path}/${out_path.baseName}_mutations.csv"), optional: true

    script:
    """
#    if [[ \$(wc -l ${assembly_path}/prot_names.txt | cut -d\$' ' -f1) == \$(wc -l ${muts_path}/stopper | cut -d\$' ' -f1) ]]; then
    if [[ \$(awk 'END {print NR}' ${assembly_path}/prot_names.txt) == \$(awk 'END {print NR}' ${muts_path}/stopper) ]]; then
        rm ${muts_path}/stopper
        for file in ${muts_path}/"${out_path.baseName}"_*_mutations.csv; do
            if [[ -f \$file ]]; then
                tail -n +2 \$file >> ${muts_path}/${out_path.baseName}_mutations.csv
                rm \$file
            fi
        done
#        sed -i 1i"SampleID;Gene;Mutation_type;Aa_change;Type_of_aa_change;Nt_mutation;Mutation_frequency" ${muts_path}/${out_path.baseName}_mutations.csv
#        sed -i '' '1i\\SampleID;Gene;Mutation_type;Aa_change;Type_of_aa_change;Nt_mutation;Mutation_frequency' ${muts_path}/${out_path.baseName}_mutations.csv
        printf "SampleID;Gene;Mutation_type;Aa_change;Type_of_aa_change;Nt_mutation;Mutation_frequency\\n\$(cat ${muts_path}/${out_path.baseName}_mutations.csv)" > temp && mv temp  ${muts_path}/${out_path.baseName}_mutations.csv
    fi
    """
}

process vizNoAnnot {
    input:
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    file ref_seq
    file muts
    val syn_muts

    output:
    stdout

    script:
    """
    python3 $projectDir/bin/viz_no_annotate.py --out-dir ${out_path} --ref-seq ${params.ref_seq} --prot-names ${assembly_path}/prot_names.txt --syn-muts ${params.syn_muts}
    """
}

process annotateAndViz {
    input:
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    file ref_seq
    file annotate
    file muts
    val syn_muts

    output:
    stdout

    script:
    """
    python3 $projectDir/bin/annotate_and_viz.py --out-dir ${out_path} --annotate ${params.annotate} --ref-seq ${params.ref_seq} --prot-names ${assembly_path}/prot_names.txt --syn-muts ${params.syn_muts}
    """
}
