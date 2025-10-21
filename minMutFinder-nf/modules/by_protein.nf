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

process protNames {
    errorStrategy 'terminate'
    input:
    path out_path
    tuple file(prot_variants_vcf), file(indelqual_prot_sam), file(samfile), file(prot_variants_vcf_gz), file(prot_variants_vcf_gz_indexed), file(vcf_file), file(prot_variants_AF_vcf_gz), file(prot_variants_AF_vcf_gz_indexed), file(prot_variants_AF_fasta)

    output:
    file "*_prot_name"

    script:
    """
    SAMPLE=${out_path.baseName}

    # Extract unique protein names from VCF file
    grep -v "#" ${vcf_file} | cut -d\$'\t' -f1 | uniq > ${out_path}/assembly/prot_names.txt

    # Check if any mutations were detected
    if [[ ! -s ${out_path}/assembly/prot_names.txt ]]; then
        echo "No mutations detected for sample \$SAMPLE"
        rm -r ${out_path}/mutations
        rm -r ${out_path}/plots
        exit 1
    fi

    # Create protein name files
    grep -v "#" ${vcf_file} | cut -d\$'\t' -f1 | uniq | awk -F "|" '{close(F); ID=\$0; gsub("^>", "", ID); F=""ID"_prot_name"} {print >> F}'
    """
}

process byProteinAnalysis {
    errorStrategy 'terminate'
    input:
    path out_path
    file protein
    val AF
    val depth
    val depth_tsv
    file samfile
    file ref_seq

    output:
    tuple path("${out_path}/mutations/stopper"), path("${out_path}/mutations/${out_path.baseName}_*_mutations.csv"), optional: true

    script:
    """
    p="${protein.baseName}"
    prot="\${p/_prot_name/}"

    echo -e "\${prot} is the protein fasta name"

    # Run the by_protein_analysis.py script
    python ${params.project_data}/bin/by_protein_analysis.py --out-dir ${out_path} --ref-seq ${ref_seq} --sample ${out_path.baseName} --prot "\${prot}" --AF ${AF} --depth ${depth} --samfile ${samfile} --depth-tsv ${depth_tsv}

    # Log the processed protein
    echo \${prot} >> ${out_path}/mutations/stopper
    """
}
