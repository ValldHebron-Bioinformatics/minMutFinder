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

process fastqProcessing {
    errorStrategy 'terminate'

    input:
    file r1
    file r2
    path out_path

    output:
    tuple file("${out_path.baseName}_qc_R1.fastq.gz"), file("${out_path.baseName}_qc_R2.fastq.gz"), file("${out_path.baseName}_qc_R1_unpaired.fastq.gz"), file("${out_path.baseName}_qc_R2_unpaired.fastq.gz")

    script:
    """
    SAMPLE=\$(basename ${out_path})

    # Filter by quality using Trimmomatic
    FQ_R1="\${SAMPLE}_qc_R1.fastq.gz"
    FQ_R2="\${SAMPLE}_qc_R2.fastq.gz"
    FQ_U_R1="\${SAMPLE}_qc_R1_unpaired.fastq.gz"
    FQ_U_R2="\${SAMPLE}_qc_R2_unpaired.fastq.gz"
    printf "\nTRIMMOMATIC: \nSoft filtering fastq files from request ID: \${SAMPLE}\n\n"
    time (trimmomatic PE -threads 64 -phred33 ${r1} ${r2} \${FQ_R1} \${FQ_U_R1} \${FQ_R2} \${FQ_U_R2} LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30)
    """
}

process mapping1_minimap2 {
    errorStrategy 'terminate'

    input:
    path out_path 
    tuple file(fqr1), file(fqr2), file(fqUPr1), file(fqUPr2)
    file ref_seq

    output:
    tuple file("${out_path.baseName}_reads-mapped_to_prot.sam"), file("${out_path.baseName}_reads-mapped_to_prot.bam")

    script:
    """
    SAMPLE=\$(basename ${out_path})
    ASSEMBLY="assembly"

    # Mapping against reference sequence using Minimap2
    echo -e "\nMINIMAP2\n->Mapping against ${ref_seq}\n"
    time (minimap2 -ax sr ${ref_seq} ${fqr1} ${fqr2} > \${SAMPLE}_reads-mapped_to_prot.sam)
    samtools view -bS \${SAMPLE}_reads-mapped_to_prot.sam | samtools sort - -o \${SAMPLE}_reads-mapped_to_prot.bam
    """
}

process mapping1_bbmap {
    errorStrategy 'terminate'

    input:
    path out_path 
    tuple file(fqr1), file(fqr2), file(fqUPr1), file(fqUPr2)
    file ref_seq

    output:
    tuple file("${out_path.baseName}_reads-mapped_to_prot.sam"), file("${out_path.baseName}_reads-mapped_to_prot.bam")

    script:
    """
    SAMPLE=\$(basename ${out_path})
    ASSEMBLY="assembly"

    # Mapping against reference sequence using Minimap2
    echo -e "\nMINIMAP2\n->Mapping against ${ref_seq}\n"
    time (bbmap.sh ref=${ref_seq} in=${fqr1} in2=${fqr2} out=\${SAMPLE}_reads-mapped_to_prot.sam sam=1.3)
    samtools view -bS \${SAMPLE}_reads-mapped_to_prot.sam | samtools sort - -o \${SAMPLE}_reads-mapped_to_prot.bam
    """
}

process variant_calling_fq {
    errorStrategy 'terminate'
    publishDir "$params.out_path/variant_calling", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.vcf"
    publishDir "$params.out_path/assembly", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.fasta"

    input:
    path out_path
    tuple file(sam), file(bam)
    val AF
    val depth
    val SB
    file ref_seq
    val threads

    output:
    tuple file("${out_path.baseName}_prot_variants.vcf"), file("${out_path.baseName}_indelqual_prot.sam"), file("samfile.sam"), file("${out_path.baseName}_prot_variants.vcf.gz"), file("${out_path.baseName}_prot_variants.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.fasta")

    script:
    """
    SAMPLE=\$(basename ${out_path})
    BAM2=${bam}

    # Call variants using LoFreq
    echo -e "\nLOFREQ\n->Variant Calling Step: Mutation detection\n"
    echo -e "lofreq indelqual --dindel -f ${ref_seq} -o "\${SAMPLE}_indelqual_prot.bam" "\$BAM2""
    lofreq indelqual --dindel -f ${ref_seq} -o "\${SAMPLE}_indelqual_prot.bam" "\$BAM2"
    INDELQUAL_BAM="\${SAMPLE}_indelqual_prot.bam"
    samtools index "\${INDELQUAL_BAM}"
    time (lofreq call-parallel --pp-threads ${threads} --call-indels -f "${ref_seq}" -o "\${SAMPLE}_prot_variants.vcf" "\${INDELQUAL_BAM}")

    samtools view -h -o "\${SAMPLE}_indelqual_prot.sam" "\${INDELQUAL_BAM}"
    awk '(\$6 != "*") && (\$12 != "") && (\$3 != "*")' "\${SAMPLE}_indelqual_prot.sam" > samfile.sam

    # Generate consensus sequence using BCFtools
    echo -e "\nBCFTOOLS\n->Variant Calling Step: Consensus sequence obtantion\n"
    bgzip -c "\${SAMPLE}_prot_variants.vcf" > "\${SAMPLE}_prot_variants.vcf.gz"
    bcftools index "\${SAMPLE}_prot_variants.vcf.gz"
    FREQ_CUTOFF=${AF}
    MIN_DEPTH=${depth}
    SB=${SB}
    bcftools filter -i'INFO/AF>='\${FREQ_CUTOFF}' && INFO/DP>='\${MIN_DEPTH}' && INFO/SB<='\${SB}'' "\${SAMPLE}_prot_variants.vcf.gz" > "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf"
    VCF_FILE_CUTOFF="\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf"
    bgzip -c "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf" > "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf.gz"
    bcftools index \${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf.gz
    cat "${ref_seq}" | bcftools consensus -s - "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf.gz" -o "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.fasta"
    """
}

process variant_calling_areads {
    errorStrategy 'terminate'
    publishDir "$params.out_path/variant_calling", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.vcf*"
    publishDir "$params.out_path/assembly", mode: 'copy', pattern: "${out_path.baseName}_prot_variants_AF${AF}.fasta"

    input:
    path out_path
    tuple file(sam), file(bam)
    val AF
    val depth
    val SB
    file ref_seq
    val threads

    output:
    tuple file("${out_path.baseName}_prot_variants.vcf"), file("${out_path.baseName}_indelqual_prot.sam"), file("samfile.sam"), file("${out_path.baseName}_prot_variants.vcf.gz"), file("${out_path.baseName}_prot_variants.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz"), file("${out_path.baseName}_prot_variants_AF${AF}.vcf.gz.csi"), file("${out_path.baseName}_prot_variants_AF${AF}.fasta")

    script:
    """
    SAMPLE=\$(basename ${out_path})

    # Call variants using LoFreq
    echo -e "\nLOFREQ\n->Variant Calling Step: Mutation detection\n"
    INDELQUAL_BAM="${bam}"
    INDELQUAL_SAM="${sam}"
    samtools index "\${INDELQUAL_BAM}"
    samtools faidx ${ref_seq}

    time (lofreq call-parallel --pp-threads ${threads} --call-indels -f "${ref_seq}" -o "\${SAMPLE}_prot_variants.vcf" "\${INDELQUAL_BAM}")

    awk '(\$6 != "*") && (\$12 != "") && (\$3 != "*")' "\${INDELQUAL_SAM}" > "samfile.sam"

    # Generate consensus sequence using BCFtools
    echo -e "\nBCFTOOLS\n->Variant Calling Step: Consensus sequence obtantion\n"
    bgzip -c "\${SAMPLE}_prot_variants.vcf" > "\${SAMPLE}_prot_variants.vcf.gz"
    bcftools index "\${SAMPLE}_prot_variants.vcf.gz"
    FREQ_CUTOFF=${AF}
    MIN_DEPTH=${depth}
    SB=${SB}
    bcftools filter -i'INFO/AF>='\${FREQ_CUTOFF}' && INFO/DP>='\${MIN_DEPTH}' && INFO/SB<='\${SB}'' "\${SAMPLE}_prot_variants.vcf.gz" > "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf"
    VCF_FILE_CUTOFF="\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf"
    bgzip -c "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf" > "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf.gz"
    bcftools index \${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf.gz
    cat ${ref_seq} | bcftools consensus -s - "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.vcf.gz" -o "\${SAMPLE}_prot_variants_AF\${FREQ_CUTOFF}.fasta"
    """
}

process mapping2_minimap2 {
    errorStrategy 'terminate'
    input:
    path out_path
    tuple file(prot_variants_vcf), file(indelqual_prot_sam), file(samfile), file(prot_variants_vcf_gz), file(prot_variants_vcf_gz_indexed), file(prot_variants_AF_vcf), file(prot_variants_AF_vcf_gz), file(prot_variants_AF_vcf_gz_indexed), file(prot_variants_AF_fasta)
    tuple file(fqr1), file(fqr2), file(fqUPr1), file(fqUPr2)

    output:
    tuple file("${out_path.baseName}_reads_mapped_consensus.sam"), file("${out_path.baseName}_reads_mapped_consensus.bam")

    script:
    """
    SAMPLE=${out_path.baseName}

    # Mapping against consensus sequence using Minimap2
    time (minimap2 -ax sr ${prot_variants_AF_fasta} ${fqr1} ${fqr2} > "\$SAMPLE"_reads_mapped_consensus.sam)
    samtools view -bS "\$SAMPLE"_reads_mapped_consensus.sam | samtools sort - -o "\$SAMPLE"_reads_mapped_consensus.bam
    """
}

process mapping2_bbmap {
    errorStrategy 'terminate'
    input:
    path out_path
    tuple file(prot_variants_vcf), file(indelqual_prot_sam), file(samfile), file(prot_variants_vcf_gz), file(prot_variants_vcf_gz_indexed), file(prot_variants_AF_vcf), file(prot_variants_AF_vcf_gz), file(prot_variants_AF_vcf_gz_indexed), file(prot_variants_AF_fasta)
    tuple file(fqr1), file(fqr2), file(fqUPr1), file(fqUPr2)

    output:
    tuple file("${out_path.baseName}_reads_mapped_consensus.sam"), file("${out_path.baseName}_reads_mapped_consensus.bam")

    script:
    """
    SAMPLE=${out_path.baseName}

    # Mapping against consensus sequence using Minimap2
    time (bbmap.sh ref=${prot_variants_AF_fasta} in=${fqr1} in2=${fqr2} out="\$SAMPLE"_reads_mapped_consensus.sam sam=1.3)
    samtools view -bS "\$SAMPLE"_reads_mapped_consensus.sam | samtools sort - -o "\$SAMPLE"_reads_mapped_consensus.bam
    """
}

process read_depth {
    errorStrategy 'terminate'
    publishDir "$params.out_path/assembly", mode: 'copy', pattern: "${out_path.baseName}_depth_consensus.tsv"
    input:
    path out_path
    tuple file(sam), file(bam)

    output:
    file "${out_path.baseName}_depth_consensus.tsv"

    script:
    """
    SAMPLE=\$(basename ${out_path})
    BAM2=${bam}

    # Calculate read depth using Samtools
    echo -e "ref\tpos\tdepth" > "\$SAMPLE"_depth_consensus.tsv
    samtools depth -a \$BAM2 >> "\$SAMPLE"_depth_consensus.tsv
    """
}

process reads_qc_metrics {
    errorStrategy 'terminate'
    publishDir "$params.out_path/qc", mode: 'copy', pattern: 'QC_metrics.csv'
    input:
    file r1
    file r2
    path out_path
    tuple file(fqr1), file(fqr2), file(fqUPr1), file(fqUPr2)
    file muts

    output:
    file "${out_path}/qc/QC_metrics.csv"

    script:
    """
    # Run the read_metrics.py script to generate QC metrics
    python3 ${params.project_data}/bin/read_metrics.py --out-dir ${out_path} --r1 ${r1} --r2 ${r2} --sample ${out_path.baseName} --r1qc ${fqr1} --r1UPqc ${fqUPr1} --r2qc ${fqr2} --r2UPqc ${fqUPr2}
    """
}
