#!/usr/bin/env python3
import sys
import os
import csv
import pandas as pd
from seq_reader import seq_reader
from grep import grep
from replacer import replacer
from variant_codon_mutation import variant_codon_mutation
from minority_mutations_finder import minority_mutations_detector
from indel_mutations_finder import indel_mutations_detector
from reads_assembler import reads_assembler

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--ref-seq':
        ref_seq = sys.argv[i + 1]
    elif sys.argv[i] == '--sample':
        SAMPLE = sys.argv[i + 1]
    elif sys.argv[i] == '--prot':
        header = sys.argv[i + 1]

FASTQ = out_dir + '/fastq'
QC_DIR = out_dir + '/qc'
qc_metrics_file = QC_DIR + '/qc_metrics.csv'
PLOTS = out_dir + '/plots'
ASSEMBLY = out_dir + '/assembly'
VARIANT_CALLING = out_dir + '/variant_calling'
MUTATIONS = out_dir + '/mutations'
FREQ_CUTOFF = str(0.05)
PROT_REF = seq_reader(ref_seq)

if os.path.isfile(ASSEMBLY + '/' + SAMPLE + '_depth_consensus.tsv'):
    df_depth_check = True
#    df_depth = pd.read_csv(ASSEMBLY + '/' + SAMPLE + '_depth_consensus.tsv', sep='\t')
else:
    df_depth_check = False
#    df_depth = False


REFERENCES = ASSEMBLY+'/references'
#if os.path.isdir(REFERENCES) is False:
#    os.mkdir(REFERENCES)

#
FREQ_CUTOFF_FASTA = seq_reader(ASSEMBLY+'/'+SAMPLE+'_prot_variants_AF-'+FREQ_CUTOFF+'.fasta')
VCF_FILE_CUTOFF = VARIANT_CALLING+'/'+SAMPLE+"_prot_variants_AF-"+FREQ_CUTOFF+".vcf"

samfile_df = open(VARIANT_CALLING + '/samfile.sam', 'r')
samfile_list = list(csv.reader(samfile_df, delimiter = "\n"))
samfile_df.close()

# df_muts_ALL = pd.DataFrame(
#     columns=[
#         "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
#         "Nt_mutation", "Mutation_frequency"
#     ]
# )

if (df_depth_check):
    depth4prot = pd.DataFrame(columns=['ref', 'pos', 'depth'])

if grep(VCF_FILE_CUTOFF, header, 'bool'):

    with open(ASSEMBLY+'/'+SAMPLE+'_prot_variants_AF-'+FREQ_CUTOFF+'_'+header+'.fasta', 'w') as fasta_seq:
        fasta_seq.write('>'+header+'\n')
        fasta_seq.write(FREQ_CUTOFF_FASTA[header]+'\n')
        fasta_seq.close()

    if (df_depth_check):
        for j in grep(ASSEMBLY + '/' + SAMPLE + '_depth_consensus.tsv', header + '\t', 'text').split('\n'):
            if len(j.split('\t')) == 3:
                depth4prot.loc[len(depth4prot)] = j.split('\t')
        depth4prot.pos, depth4prot.depth = pd.to_numeric(depth4prot.pos), pd.to_numeric(depth4prot.depth)

        depth4prot_corrected = pd.DataFrame(columns=['ref', 'pos', 'depth'])
        for p in range(1, len(FREQ_CUTOFF_FASTA[header])+1):
            if depth4prot[depth4prot.pos == p].empty:
                row = header + ' ' + str(p) + ' ' + '0'
            else:
                row = depth4prot[depth4prot.pos == p].to_string(index=False, header=False)
            depth4prot_corrected.loc[len(depth4prot_corrected)] = row.split(' ')

        depth4prot_corrected.pos = pd.to_numeric(depth4prot_corrected.pos)
        depth4prot_corrected.depth = pd.to_numeric(depth4prot_corrected.depth)

        depth4prot_corrected.to_csv(ASSEMBLY+'/'+header+'_depth.tsv', sep='\t', index=False)
        file = pd.read_csv(ASSEMBLY+'/'+header+'_depth.tsv', sep="\t")


        sequence_file = open(ASSEMBLY+'/'+SAMPLE+'_prot_variants_AF-'+FREQ_CUTOFF+'_'+header+'.fasta', "r")

        # Store de fasta in one line sequence
        sequence = [i.rstrip() for i in sequence_file.readlines()]
        tag = sequence.pop(0)
        sequence = "".join(sequence)
        sequence_file.close()

        for index, row in file.iterrows():
            # N postions replacer
            #print('\n')
            #print('Replacing positions with coverage below 20 with "N" in consensus sequence...')
            #print('\n')
            if int(row['depth']) < 20:
                sequence = replacer(sequence, "N", int(row['pos'])-1)

            else:
                continue

            sequence_format = [sequence[n:n + 70] for n in range(0, len(sequence), 70)]
            outfile = open(ASSEMBLY + '/' + SAMPLE + '_prot_variants_AF-' + FREQ_CUTOFF + '_'+header + '.fasta', 'w')

            outfile.write(tag)
            outfile.write('\n')
            for seq in sequence_format:
                outfile.write(seq)
                outfile.write('\n')
        ##
    FREQ_CUTOFF_FASTA_prot = seq_reader(ASSEMBLY+'/'+SAMPLE+'_prot_variants_AF-'+FREQ_CUTOFF+'_'+header+'.fasta')

    # SEQ METRICS
    length_assembly = len(FREQ_CUTOFF_FASTA_prot[header])

    if (df_depth_check):
        n_percentage = ((FREQ_CUTOFF_FASTA_prot[header].count('N')) / length_assembly)*100
        n_percentage_clean = str(round(n_percentage, 2))+'%'

    ref_seq, cons_seq = PROT_REF[header].upper(), FREQ_CUTOFF_FASTA_prot[header].upper()

    cons_seq = cons_seq.upper()
    cons_seq = cons_seq.replace('N', '')
    coverage = (len(cons_seq) / len(ref_seq)) * 100
    coverage_clean = str(round(coverage, 2)) + '%'
    with open(qc_metrics_file, 'a') as infile:
        infile.write(SAMPLE+';'+'length_'+header+';'+str(length_assembly)+'\n')
        if (df_depth_check):
            infile.write(SAMPLE+';'+'n_percentage_'+header+';'+str(n_percentage_clean)+'\n')
        infile.write(SAMPLE+';'+'coverage_'+header+';'+str(coverage_clean)+'\n')
        if (df_depth_check):
            infile.write(SAMPLE+';'+'median_'+header+';'+str(depth4prot_corrected.depth.quantile([0.5]).iloc[0])+'\n')
            infile.write(SAMPLE+';'+'Q1_'+header+';'+str(depth4prot_corrected.depth.quantile([0.25]).iloc[0])+'\n')
            infile.write(SAMPLE+';'+'Q3_'+header+';'+str(depth4prot_corrected.depth.quantile([0.75]).iloc[0])+'\n')
    infile.close()

    qc_metrics_file_df = pd.read_csv(qc_metrics_file, sep=';', header=None)
    qc_metrics_protein_file = pd.DataFrame(columns=["sample", "test", "score"])

    for q in grep(qc_metrics_file_df, header, 'text'):
        q_df = pd.DataFrame(
            grep(
                qc_metrics_file_df, header, 'text'
            )
        )[
            pd.DataFrame(grep(qc_metrics_file_df, header, 'text')).index == q
        ].to_string(index=False, header=None).split(' ')

    if len(q_df) == len(qc_metrics_protein_file.columns):
        qc_metrics_protein_file.loc[len(qc_metrics_protein_file)] = q_df
    qc_metrics_protein_file.to_csv(QC_DIR+'/'+header+'_qc_metrics.csv', sep=';', index=False)

    print('\nBeggining analyisis for protein '+header+'\n')
    vcf_header = grep(VCF_FILE_CUTOFF, '#C', 'text').replace('\n', '').replace('#', '').split('\t')

    vcf_prot = pd.DataFrame(columns=vcf_header)
    for v in grep(VCF_FILE_CUTOFF, header, 'text').split('\n'):
        if len(v.split('\t')) == len(vcf_header):
            vcf_prot.loc[len(vcf_prot)] = v.split('\t')

    samfile_prot, max_read_len = [], 0
    for list_line in samfile_list:
        line = list_line[0]
        if header == line.split('\t')[2]:
            if max_read_len < len(line.split('\t')[9]):
                max_read_len = len(line.split('\t')[9])
            samfile_prot.append(line)
    max_read_aa_len = max_read_len//3
    if (max_read_len % 3) != 0:
        max_read_aa_len += 1

    df_sam = pd.DataFrame(pd.DataFrame(samfile_prot)[0].str.split(pat='\t').to_list())
    df_sam.to_csv(VARIANT_CALLING+'/samfile_'+header+'.tsv', sep='\t', index=False, header=None)

    df_sam[3] = pd.to_numeric(df_sam[3])
    VCF_FILE = ASSEMBLY+'/'+SAMPLE+'_'+header+'_prot_variants_AF-'+FREQ_CUTOFF+'.vcf'
    SAM_TO_LOOP = VARIANT_CALLING+'/sam_to_loop_'+header+'.sam'
    READS_MULTIFASTA = VARIANT_CALLING+'/reads_multifasta_'+header+'.fasta'

    if vcf_prot.size > 0:
        vcf_prot.POS = pd.to_numeric(vcf_prot.POS)
        vcf_prot['AA_POS'] = ((vcf_prot.POS - 1) // 3) + 1
        vcf_prot['AA_INDEX'] = ((vcf_prot.POS - 1) % 3)
        AF = pd.DataFrame(
            pd.DataFrame(
                pd.DataFrame(
                    pd.DataFrame(
                        vcf_prot.INFO.str.split(pat=';')
                    )['INFO'].to_list()
                )[1].str.split(pat='=')
            )[1].to_list()
        )[1]

        vcf_prot['AF'] = AF
        vcf_prot.AF = pd.to_numeric(vcf_prot.AF)
        
        indels = vcf_prot[(vcf_prot.ALT.str.len() > 1) | (vcf_prot.REF.str.len() > 1)]
        alone_aa =  vcf_prot[(vcf_prot.REF.str.len() == 1) & (vcf_prot.ALT.str.len() == 1)]
        # frameshift = vcf_prot[vcf_prot.ALT.str.len() > 1]

        alone_aa = alone_aa.drop_duplicates(subset=['AA_POS'], keep=False)
        alone_aa = alone_aa.reset_index(drop=True)  # ready for analysis

        vcf_prot_for_reps = pd.concat([vcf_prot, alone_aa]).drop_duplicates(keep=False).reset_index(drop=True)

        repeated_aa = vcf_prot_for_reps[vcf_prot_for_reps.duplicated(subset='AA_POS', keep=False)]
        repeated_aa = repeated_aa.reset_index(drop=True)

        fixed_muts = pd.DataFrame(columns=repeated_aa.columns)
        non_fixed_muts = pd.DataFrame(columns=repeated_aa.columns)

        for rep in repeated_aa.AA_POS.unique():
            AFs, isit = repeated_aa[repeated_aa.AA_POS == rep].AF, False
            for n in AFs:
                if n >= 0.99:
                    isit = True
            if isit:
                fixed_mut = repeated_aa[repeated_aa.AA_POS == rep]
                if fixed_muts.size > 0:
                    fixed_muts = pd.concat([fixed_muts, fixed_mut])
                else:
                    fixed_muts = fixed_mut
                print(str(rep) + ' is a fixed mutation')

            else:
                non_fixed_mut = repeated_aa[repeated_aa.AA_POS == rep]
                if non_fixed_muts.size > 0:
                    non_fixed_muts = pd.concat([non_fixed_muts, non_fixed_mut])
                else:
                    non_fixed_muts = non_fixed_mut
                print(str(rep) + ' is not a fixed mutation')

        fixed_muts = fixed_muts.reset_index(drop=True)  # ready for analysis
        non_fixed_muts = non_fixed_muts.reset_index(drop=True)  # ready for analysis
        if alone_aa.size > 0:
            muts = []
            for n_aa in alone_aa.AA_POS.unique():
                muts_aux, n_aa_df = [], alone_aa[alone_aa.AA_POS == n_aa]
                n_aa_df = n_aa_df.reset_index(drop=True)
                muts_aux = minority_mutations_detector(n_aa_df, ref_seq, SAMPLE, header)
                muts.append(muts_aux)
            print(muts)
            dfmuts = pd.DataFrame(
                muts, columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                    "Nt_mutation", "Mutation_frequency"
                ]
            )

            dfmuts.to_csv(
                MUTATIONS + '/' + SAMPLE + '_' + header + '_alone_mutations.csv', sep=';', index=False
            )

        if indels.size > 0:
            muts = pd.DataFrame(columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                    "Nt_mutation", "Mutation_frequency"
                ])
            for n_aa in indels.AA_POS.unique():
                muts_aux, n_aa_df = [], indels[indels.AA_POS == n_aa]
                n_aa_df = n_aa_df.reset_index(drop=True)
                muts_aux = indel_mutations_detector(n_aa_df, ref_seq, SAMPLE, header)
                muts_aux = muts_aux.reset_index(drop=True)
                muts = pd.concat([muts, muts_aux]).reset_index(drop=True)
            print(muts)

            muts.to_csv(
                MUTATIONS + '/' + SAMPLE + '_' + header + '_indel_mutations.csv', sep=';', index=False
            )

        if fixed_muts.size > 0:
            muts = pd.DataFrame(columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                    "Nt_mutation", "Mutation_frequency"
                ])
            for n_aa in fixed_muts.AA_POS.unique():
                muts_aux, n_aa_df = [], fixed_muts[fixed_muts.AA_POS == n_aa]
                n_aa_df = n_aa_df.reset_index(drop=True)
                muts_aux = variant_codon_mutation(n_aa_df, ref_seq, "Fixed", n_aa, SAMPLE, header)
                muts_aux = muts_aux.reset_index(drop=True)
                muts = pd.concat([muts, muts_aux]).reset_index(drop=True)

            dfmuts = pd.DataFrame(
                muts, columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                    "Nt_mutation", "Mutation_frequency"
                ]
            )

            dfmuts.to_csv(MUTATIONS + '/' + SAMPLE + '_' + header + '_fixed_mutations.csv', sep=';', index=False)

        if non_fixed_muts.size > 0:
            min_aa, max_aa = min(non_fixed_muts.AA_POS), max(non_fixed_muts.AA_POS)
            minimum, maximum = min_aa * 3, max_aa * 3

            print('Maximum read length for protein ' + header + ' is -> ' + str(max_read_len))
            if (maximum - minimum) > max_read_len:
                current = ''
                for aa_ranges in non_fixed_muts.AA_POS.unique():
                    previous = current
                    current, minimum = aa_ranges, aa_ranges * 3
                    if previous != '':
                        range_max = minimum + 3
                        if (current - previous) >= max_read_aa_len or (current - previous) <= - max_read_aa_len:
                            range_min = minimum - max_read_len
                        else:
                            nt_previous = previous * 3
                            range_min = nt_previous
                            # take from samfile regioninside rangemax and rangemin

                    else:
                        range_min = minimum - max_read_len
                        range_max = minimum + max_read_len
                        # take from samfile regioninside rangemax and rangemin
    
                    sam_to_loop = df_sam[df_sam[3] >= range_min][df_sam[df_sam[3] >= range_min][3] <= range_max]
                    sam_to_loop = sam_to_loop.reset_index(drop=True)
                    sam_to_loop.to_csv(SAM_TO_LOOP, index=False, sep='\t', header=None)
                    minimum += 203

                    reads_fasta = reads_assembler(SAM_TO_LOOP, '')
                    with open(READS_MULTIFASTA, "a+") as ofile:
                        ofile.write(reads_fasta + '\n')
                    # take from samfile regioninside rangemax and rangemin
            else:
                range_max, range_min = maximum + 3, minimum - max_read_len
                print(range_min)
                print(range_max)
                sam_to_loop = df_sam[df_sam[3] >= range_min][df_sam[df_sam[3] >= range_min][3] <= range_max]
                sam_to_loop = sam_to_loop.reset_index(drop=True)
                sam_to_loop.to_csv(SAM_TO_LOOP, index=False, sep='\t', header=None)
                reads_fasta = reads_assembler(SAM_TO_LOOP, '')
                with open(READS_MULTIFASTA, "a+") as ofile:
                    ofile.write(reads_fasta+'\n')
                # take from samfile regioninside rangemax and rangemin
                # ASSEMBLY_PY
            reads = open(READS_MULTIFASTA, 'r')
            reads_df = pd.DataFrame(list(csv.reader(reads, delimiter="\n")))
            reads.close()

            just_reads = []
            for i in range(0, reads_df.shape[0] - 1):
                if reads_df[0][i] is not None:
                    if '>' not in reads_df[0][i]:
                        just_reads.append(reads_df.iloc[[i]].to_string(index=False, header=None))
                    else:
                        continue
                else:
                    continue

            just_reads_df = pd.DataFrame(just_reads)

            # Extract strings from the column and split them into groups of 3 characters
            just_reads_df_sep = just_reads_df[0].str.extractall('(.{1,3})')[0].unstack().fillna('---')

            muts = pd.DataFrame(columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                    "Nt_mutation", "Mutation_frequency"
                ])
            for n_aa in non_fixed_muts.AA_POS.unique():
                muts_aux, n_aa_df = [], non_fixed_muts[non_fixed_muts.AA_POS == n_aa]
                n_aa_df = n_aa_df.reset_index(drop=True)
                muts_aux = variant_codon_mutation(n_aa_df, ref_seq, just_reads_df_sep, n_aa, SAMPLE, header)
                muts_aux = muts_aux.reset_index(drop=True)
                muts = pd.concat([muts, muts_aux]).reset_index(drop=True)

            dfmuts = pd.DataFrame(
                muts, columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change",
                    "Type_of_aa_change ", "Nt_mutation", "Mutation_frequency"
                ]
            )

            # if df_muts_ALL.size > 0:
            #     df_muts_ALL = pd.concat([df_muts_ALL, dfmuts])
            # else:
            #     df_muts_ALL = dfmuts
            dfmuts.to_csv(MUTATIONS + '/' + SAMPLE + '_' + header + '_not_fixed_mutations.csv', sep=';', index=False)

#df_muts_ALL.to_csv(MUTATIONS + '/' + SAMPLE + '_mutations.csv', sep=';', index=False, header=None)
