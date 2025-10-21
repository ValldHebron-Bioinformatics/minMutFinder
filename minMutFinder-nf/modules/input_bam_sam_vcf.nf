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

process inVcf {
    errorStrategy 'terminate'
    publishDir "$params.out_path/variant_calling", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.vcf*"
    publishDir "$params.out_path/assembly", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.fasta"

    input:
    file vcf
    path out_path
    tuple file(sam), file(bam)
    val AF
    val depth
    val SB
    file ref_seq

    output:
    tuple path("${out_path}/variant_calling/${out_path.baseName}_prot_variants.vcf"), file("${out_path.baseName}_indelqual_prot.sam"), file("samfile.sam"), file("${out_path.baseName}_prot_variants.vcf.gz"), file("${out_path.baseName}_prot_variants.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.fasta")

    script:
    """
    SAMPLE=${out_path.baseName}

    # Filter SAM file
    awk '(\$6 != "*") && (\$12 != "") && (\$3 != "*")' ${sam} > "samfile.sam"

    # Consensus sequence
    echo -e "\nBCFTOOLS\n->Variant Calling Step: Consensus sequence obtantion\n"
    if [[ ${params.vcf} == *".vcf"* ]]; then
        python3 ${params.project_data}/bin/in_vcf.py --out-dir ${params.out_path} --vcf ${params.vcf}
    else
        python3 ${params.project_data}/bin/tsv2vcf.py --out-dir ${params.out_path} --tsv ${params.vcf}
    fi

    # Compress and index VCF file
    bgzip -c ${vcf} > ""\${SAMPLE}"_prot_variants.vcf.gz"
    bcftools index ""\${SAMPLE}"_prot_variants.vcf.gz"

    # Filter VCF file by allele frequency and depth
    FREQ_CUTOFF=${AF}
    MIN_DEPTH=${depth}
    SB=${SB}
    bcftools filter -i'INFO/AF>='\${FREQ_CUTOFF}' && INFO/DP>='\${MIN_DEPTH}' && INFO/SB<='\${SB}'' ""\${SAMPLE}"_prot_variants.vcf.gz" > ""\${SAMPLE}"_prot_variants_AF"\${FREQ_CUTOFF}".vcf"
    VCF_FILE_CUTOFF=""\${SAMPLE}"_prot_variants_AF\${FREQ_CUTOFF}.vcf"
    bgzip -c "\${VCF_FILE_CUTOFF}" > "\${VCF_FILE_CUTOFF}.gz"
    bcftools index "\${VCF_FILE_CUTOFF}.gz"

    # Generate consensus sequence
    cat ${ref_seq} | bcftools consensus -s - "\${VCF_FILE_CUTOFF}".gz -o ""\${SAMPLE}"_prot_variants_AF"\${FREQ_CUTOFF}".fasta"
    """
}

process inVcf_noSB {
    errorStrategy 'terminate'
    publishDir "$params.out_path/variant_calling", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.vcf*"
    publishDir "$params.out_path/assembly", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.fasta"

    input:
    file vcf
    path out_path
    tuple file(sam), file(bam)
    val AF
    val depth
    file ref_seq

    output:
    tuple path("${out_path}/variant_calling/${out_path.baseName}_prot_variants.vcf"), file("${out_path.baseName}_indelqual_prot.sam"), file("samfile.sam"), file("${out_path.baseName}_prot_variants.vcf.gz"), file("${out_path.baseName}_prot_variants.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.fasta")

    script:
    """
    SAMPLE=${out_path.baseName}

    # Filter SAM file
    awk '(\$6 != "*") && (\$12 != "") && (\$3 != "*")' ${sam} > "samfile.sam"

    # Consensus sequence
    echo -e "\nBCFTOOLS\n->Variant Calling Step: Consensus sequence obtantion\n"
    if [[ ${vcf} == *".vcf"* ]]; then
        python3 ${params.project_data}/bin/in_vcf.py --out-dir ${out_path} --vcf ${vcf}
    else
        python3 ${params.project_data}/bin/tsv2vcf.py --out-dir ${out_path} --tsv ${vcf}
    fi

    # Compress and index VCF file
    bgzip -c "${out_path}/variant_calling/${out_path.baseName}_prot_variants.vcf" > ""\${SAMPLE}"_prot_variants.vcf.gz"
    bcftools index ""\${SAMPLE}"_prot_variants.vcf.gz"

    # Filter VCF file by allele frequency and depth
    FREQ_CUTOFF=${AF}
    MIN_DEPTH=${depth}
    bcftools filter -i'INFO/AF>='\${FREQ_CUTOFF}' && INFO/DP>='\${MIN_DEPTH}'' ""\${SAMPLE}"_prot_variants.vcf.gz" > ""\${SAMPLE}"_prot_variants_AF"\${FREQ_CUTOFF}".vcf"
    VCF_FILE_CUTOFF=""\${SAMPLE}"_prot_variants_AF\${FREQ_CUTOFF}.vcf"
    bgzip -c "\${VCF_FILE_CUTOFF}" > "\${VCF_FILE_CUTOFF}.gz"
    bcftools index "\${VCF_FILE_CUTOFF}.gz"

    # Generate consensus sequence
    cat ${ref_seq} | bcftools consensus -s - "\${VCF_FILE_CUTOFF}".gz -o ""\${SAMPLE}"_prot_variants_AF"\${FREQ_CUTOFF}".fasta"
    """
}

process inAlignedReads {
    errorStrategy 'terminate'
    input:
    file areads
    path out_path

    output:
    tuple path("${out_path}/variant_calling/${out_path.baseName}_indelqual_prot.sam"), path("${out_path}/variant_calling/${out_path.baseName}_indelqual_prot.bam")

    script:
    """
    python3 ${params.project_data}/bin/sam_bam_processing.py --out-dir ${out_path} --areads ${areads}
    """
}
