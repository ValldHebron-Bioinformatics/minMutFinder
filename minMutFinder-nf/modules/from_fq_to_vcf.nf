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

#!/usr/bin/env nextflow

process fastqProcessing {
    input:
    file r1
    file r2
    file ref_seq
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path

    output:
    file "${assembly_path}/ref_seq.fasta"
    file "${assembly_path}/ref_seq.fasta.fai"
    file "${fq_path}/${out_path.baseName}_qc_R1.fastq.gz"
    file "${fq_path}/${out_path.baseName}_qc_R2.fastq.gz"
    file "${fq_path}/${out_path.baseName}_qc_R1_unpaired.fastq.gz"
    file "${fq_path}/${out_path.baseName}_qc_R2_unpaired.fastq.gz"

    script:
    """
    DIR_SAMPLE=${out_path}
    FASTQ=${fq_path}
    ASSEMBLY=${assembly_path}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    echo "\${SAMPLE}"

    R1=${params.r1}
    R2=${params.r2}
    cp ${params.ref_seq} "\${ASSEMBLY}/ref_seq.fasta"
    samtools faidx "\${ASSEMBLY}/ref_seq.fasta"
    # Filter by quality
    # Filter by Trimmomatic
    FQ_R1="\${FASTQ}"/"\${SAMPLE}"_qc_R1.fastq.gz; FQ_R2="\${FASTQ}"/"\${SAMPLE}"_qc_R2.fastq.gz
    printf "\nTRIMMOMATIC: \nSoft filtering fastq files from request ID: \${SAMPLE}\n\n"
    time (trimmomatic PE -threads 64 -phred33 "\$R1" "\$R2" "\${FQ_R1}" "\${FASTQ}"/"\${SAMPLE}"_qc_R1_unpaired.fastq.gz "\${FQ_R2}" "\${FASTQ}"/"\${SAMPLE}"_qc_R2_unpaired.fastq.gz LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30)
    """
}

process mapping1 {
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
    file ref_seq_fai
    file fqr1
    file fqr2
    file fqUPr1
    file fqUPr2

    output:
    file "${assembly_path}/${out_path.baseName}_reads-mapped_to_prot.sam"
    file "${assembly_path}/${out_path.baseName}_reads-mapped_to_prot.bam"

    script:
    """
    DIR_SAMPLE=${out_path}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${assembly_path}
    FQ_R1=${fqr1}
    FQ_R2=${fqr2}
    echo -e "\nMINIMAP2\n->Mapping against "\${ASSEMBLY}/ref_seq.fasta"\n"
    time (minimap2 -ax sr "\${ASSEMBLY}/ref_seq.fasta" \$FQ_R1 \$FQ_R2 > "\${ASSEMBLY}"/"\${SAMPLE}"_reads-mapped_to_prot.sam)
    samtools view -bS "\${ASSEMBLY}"/"\${SAMPLE}"_reads-mapped_to_prot.sam | samtools sort - -o "\${ASSEMBLY}"/"\${SAMPLE}"_reads-mapped_to_prot.bam
    """

}

process variant_calling {
    input:
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    file sam
    file bam


    output:
    file "${vcf_path}/${out_path.baseName}_indelqual_prot.bam"
    file "${vcf_path}/${out_path.baseName}_indelqual_prot.bam.bai"
    file "${vcf_path}/${out_path.baseName}_prot_variants.vcf"
    file "${vcf_path}/${out_path.baseName}_indelqual_prot.sam"
    file "${vcf_path}/samfile.sam"
    file "${vcf_path}/${out_path.baseName}_prot_variants.vcf.gz"
    file "${vcf_path}/${out_path.baseName}_prot_variants.vcf.gz.csi"
    file "${vcf_path}/${out_path.baseName}_prot_variants_AF-0.05.vcf"
    file "${vcf_path}/${out_path.baseName}_prot_variants_AF-0.05.vcf.gz"
    file "${vcf_path}/${out_path.baseName}_prot_variants_AF-0.05.vcf.gz.csi"
    file "${assembly_path}/${out_path.baseName}_prot_variants_AF-0.05.fasta"
    
    script:
    """
    DIR_SAMPLE=${out_path}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${assembly_path}
    VARIANT_CALLING=${vcf_path}
    BAM2=${bam}

    # Call variants
    # tomake lofreq work 
    echo -e "\nLOFREQ\n->Variant Calling Step: Mutation detection\n"
    
    lofreq indelqual --dindel -f "\${ASSEMBLY}/ref_seq.fasta" -o "\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.bam "\$BAM2"
    INDELQUAL_BAM="\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.bam
    samtools index "\${INDELQUAL_BAM}"
    # tomake lofreq work
    time (lofreq call-parallel --pp-threads 64 --call-indels -f "\${ASSEMBLY}/ref_seq.fasta" -o "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf "\${INDELQUAL_BAM}")
    
    samtools index "\${INDELQUAL_BAM}"
    samtools view -h -o "\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.sam "\${INDELQUAL_BAM}"
   
    awk '(\$6 != "*") && (\$12 != "") && (\$3 != "*")' "\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.sam > "\${VARIANT_CALLING}"/samfile.sam
   
    # Consensus sequence
    echo -e "\nBCFTOOLS\n->Variant Calling Step: Consensus sequence obtantion\n"
   
    bgzip -c "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf > "\${VARIANT_CALLING}"/"\$SAMPLE"_prot_variants.vcf.gz
    bcftools index "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf.gz
    FREQ_CUTOFF=0.05
    bcftools filter -i'AF>='\${FREQ_CUTOFF}' & DP>=20' "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf.gz > "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf
    VCF_FILE_CUTOFF="\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf
    bgzip -c "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf > "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf.gz
    bcftools index \${VARIANT_CALLING}/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf.gz
    cat "\${ASSEMBLY}/ref_seq.fasta" | bcftools consensus -s - "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf.gz -o "\${ASSEMBLY}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".fasta
    # bcftools consensus -f "\${ASSEMBLY}/ref_seq.fasta" "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf.gz -o "\${ASSEMBLY}"/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".fasta
    """

}

process mapping2 {
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
    file prot_variants_AF05_vcf
    file prot_variants_AF05_vcf_gz
    file prot_variants_AF05_vcf_gz_indexed
    file prot_variants_AF05_fasta

    output:
    file "${assembly_path}/${out_path.baseName}_reads_mapped_consensus.sam"
    file "${assembly_path}/${out_path.baseName}_reads_mapped_consensus.bam"
    
    script:
    """
    DIR_SAMPLE=${out_path}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${assembly_path}
    VARIANT_CALLING=${vcf_path}
    FASTQ=${fq_path}
    FREQ_CUTOFF=0.05
    FQ_R1="\${FASTQ}"/"\${SAMPLE}"_qc_R1.fastq.gz; FQ_R2="\${FASTQ}"/"\${SAMPLE}"_qc_R2.fastq.gz
    time (minimap2 -ax sr "${prot_variants_AF05_fasta}" \$FQ_R1 \$FQ_R2 > \$ASSEMBLY/"\$SAMPLE"_reads_mapped_consensus.sam)
    samtools view -bS \$ASSEMBLY/"\$SAMPLE"_reads_mapped_consensus.sam | samtools sort - -o \$ASSEMBLY/"\$SAMPLE"_reads_mapped_consensus.bam
    """


}


process read_depth {
    input:
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    file sam
    file bam


    output:
    file "${assembly_path}/${out_path.baseName}_depth_consensus.tsv"
    
    script:
    """
    DIR_SAMPLE=${out_path}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${assembly_path}
    VARIANT_CALLING=${vcf_path}
    FREQ_CUTOFF=0.05
    BAM2=${bam}
    echo -e "ref\tpos\tdepth" > \$ASSEMBLY/"\$SAMPLE"_depth_consensus.tsv
    samtools depth -a \$BAM2 >> \$ASSEMBLY/"\$SAMPLE"_depth_consensus.tsv
    """

}

process reads_qc_plot {
    input:
    file r1
    file r2
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    file ref_seq
    file ref_seq_fai
    file fqr1
    file fqr2
    file fqUPr1
    file fqUPr2

    output:
    file "${qc_path}/qc_metrics.csv"
    file "${plots_path}/qc_metrics.svg"

    script:
    """
    python3 $projectDir/bin/read_metrics.py --out-dir ${out_path} --r1 ${params.r1} --r2 ${params.r2} --sample ${out_path.baseName} --r1qc ${fqr1} --r1UPqc ${fqUPr1} --r2qc ${fqr2} --r2UPqc ${fqUPr2}
    """
}

