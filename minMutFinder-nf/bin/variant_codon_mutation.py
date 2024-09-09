#!/usr/bin/env python3
import pandas as pd


from replacer_vcm import replacer_vcm
from all_muts_finder import all_muts_finder

aa_classes = {
    'A': 'Hydrophobic', 'V': 'Hydrophobic',
    'L': 'Hydrophobic', 'M': 'Hydrophobic',
    'I': 'Hydrophobic', 'S': 'Polar_non_charged',
    'T': 'Polar_non_charged', 'N': 'Polar_non_charged',
    'Q': 'Polar_non_charged', 'G': 'Special_case',
    'C': 'Special_case', 'P': 'Special_case',
    'U': 'Special_case', 'F': 'Aromatic_Hydrophobic',
    'Y': 'Aromatic_Hydrophobic', 'W': 'Aromatic_Hydrophobic',
    'K': 'Positively_charged', 'R': 'Positively_charged',
    'H': 'Positively_charged', 'D': 'Positively_charged',
    'E': 'Positively_charged', '*': 'Stop'
}

codon_code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def variant_codon_mutation(codon_data, ref_seq, reads, aa_pos, sample_id, gene):
    muts = []

    # Data loading
    codon, fasta, aa_pos = codon_data, ref_seq, int(aa_pos)

    reference_codon = fasta[(aa_pos * 3) - 3:(aa_pos * 3)]
    reference_codon = reference_codon.upper()

    if type(reads) is pd.DataFrame:
        r, codon_dic = [0, 1, 2], {}
        for i in r:
            codon_dic[i] = codon["ALT"][codon["AA_INDEX"] == i].to_string(index=False).split('\n')

        muts_list = []
        for i in r:
            new_r = [0, 1, 2]
            new_r.pop(i)
            for j in codon_dic[i]:
                mutated_codon = replacer_vcm(reference_codon, j, i)
                if "Series" in mutated_codon:
                    continue
                else:
                    if mutated_codon not in muts_list:
                        muts_list.append(mutated_codon)
                    all_muts_finder(mutated_codon, new_r, codon_dic, muts_list)

        reads_seq = reads
        count, col = 0, reads_seq[aa_pos - 1][reads_seq[aa_pos - 1].notna()]
        for i in col:
            if "-" not in i:
                count += 1

        mut_dic, mut_aa_dic = {}, {}
        for mut in muts_list:
            if reference_codon != mut:
                n_reads_mut = reads_seq[reads_seq[aa_pos-1] == mut].shape[0]
                percentage, gene_pos = (n_reads_mut / count), gene  # *100
                if percentage >= 0.05:

                    gene_pos, mut_dic[mut], prot_sub = gene, percentage, codon_code[mut]
                    prot_ref, nt_mutation = codon_code[reference_codon], ''

                    mut_aa_dic[prot_sub] = prot_ref
                    for i in r:
                        if reference_codon[i] != mut[i]:
                            if nt_mutation == '':
                                nt_mutation += str(str(reference_codon[i]) + str(i + (aa_pos * 3) - 2) + str(mut[i]))
                            else:
                                nt_mutation += str(", " + str(reference_codon[i]) + str(i + (aa_pos * 3) - 2) + str(mut[i]))

                    if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                        aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                    elif (aa_classes[prot_ref] == "Special_case") and (aa_classes[prot_sub] == "Special_case"):
                        if (prot_ref != prot_sub):
                            aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                        else:
                            aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                    else:
                        aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                    if prot_ref != prot_sub:
                        mut_aux = [
                            sample_id, gene_pos, "MUTATION",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage)
                        ]

                    else:
                        mut_aux = [
                            sample_id, gene_pos, "NO_MUTATION",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage)
                        ]
                    muts.append(mut_aux)

    else:
        prot_ref = codon_code[reference_codon]
        r = [0, 1, 2]
        codon_dic = {}
        for i in r:
            codon_dic[i] = codon["ALT"][codon["AA_INDEX"] == i].to_string(index=False).split('\n')

        fixed_nts = reference_codon
        r_fixed = []
        r_not_fixed = []
        for index, row in codon.iterrows():
            ALT = codon.ALT.iloc[index]
            AA_INDEX = codon.AA_INDEX.iloc[index]
            if float(codon.AF.iloc[index]) >= 0.99:
                fixed_nts = replacer_vcm(fixed_nts, str(ALT), int(AA_INDEX))
                r_fixed.append(AA_INDEX)
                fixed_mut_freq = float(codon.AF.iloc[index])
            else:
                r_not_fixed.append(AA_INDEX)
        gene_pos = gene

        muts_list = []
        for w in r:
            new_r = [0, 1, 2]
            new_r.pop(w)
            for i in r_not_fixed:
                for j in codon_dic[i]:
                    mutated_codon = replacer_vcm(fixed_nts, j, i)
                    if "Series" in mutated_codon:
                        continue
                    else:
                        if mutated_codon not in muts_list:
                            muts_list.append(mutated_codon)
                            a = codon[codon["ALT"] == j]
                            b = a[a["AA_INDEX"] == i]
                            freq_mut = b["AF"].to_string(index=False).split('\n')[0]
                            freq_non_fixed = float(freq_mut)
                            fixed_mut_freq = fixed_mut_freq - freq_non_fixed
                            prot_sub = codon_code[mutated_codon]

                            if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                                aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                            elif (aa_classes[prot_ref] == "Special_case") and (aa_classes[prot_sub] == "Special_case"):
                                if (prot_ref != prot_sub):
                                    aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]
                                else:
                                    aa_change = "No_Change, stayed " + aa_classes[prot_ref]
                            else:
                                aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                            nt_mutation = ''
                            for y in r:
                                if reference_codon[y] != mutated_codon[y]:
                                    if nt_mutation == '':
                                        nt_mutation += str(
                                            str(reference_codon[y]) + str(y + (aa_pos * 3) - 2) + str(mutated_codon[y])
                                        )
                                    else:
                                        nt_mutation += str(
                                            ", " + str(
                                                reference_codon[y]
                                            ) + str(
                                                y + (aa_pos * 3) - 2
                                            ) + str(
                                                mutated_codon[y]
                                            )
                                        )

                            percentage = freq_non_fixed
                            if percentage >= 0.05:
                                if prot_ref != prot_sub:
                                    mut_aux = [
                                        sample_id, gene_pos, "MUTATION",
                                        prot_ref + str(aa_pos) + prot_sub,
                                        aa_change, nt_mutation, str(percentage)
                                    ]
                                else:
                                    mut_aux = [
                                        sample_id, gene_pos, "NO_MUTATION",
                                        prot_ref + str(aa_pos) + prot_sub,
                                        aa_change, nt_mutation, str(percentage)
                                    ]
                                muts.append(mut_aux)
                else:
                    continue

        nt_fixed_list = []
        for i in r_fixed:
            prot_sub = codon_code[fixed_nts]
            if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]
            elif (aa_classes[prot_ref] == "Special_case") and (aa_classes[prot_sub] == "Special_case"):
                if (prot_ref != prot_sub):
                    aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]
                else:
                    aa_change = "No_Change, stayed " + aa_classes[prot_ref]
            else:
                aa_change = "No_Change, stayed " + aa_classes[prot_ref]

            nt_mutation = ''
            for i in r:
                if reference_codon[i] != fixed_nts[i]:
                    if nt_mutation == '':
                        nt_mutation += str(
                            str(
                                reference_codon[i]
                            ) + str(
                                i + (aa_pos * 3) - 2
                            ) + str(
                                fixed_nts[i]
                            )
                        )
                    else:
                        nt_mutation += str(
                            ", " + str(
                                reference_codon[i]
                            ) + str(
                                i + (aa_pos * 3) - 2
                            ) + str(
                                fixed_nts[i]
                            )
                        )

            if nt_mutation not in nt_fixed_list:
                nt_fixed_list.append(nt_mutation)
                percentage = fixed_mut_freq

                if percentage >= 0.05:
                    if prot_ref != prot_sub:
                        mut_aux = [
                            sample_id, gene_pos, "MUTATION",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage)
                        ]
                    else:
                        mut_aux = [
                            sample_id, gene_pos, "NO_MUTATION",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage)
                        ]
                    muts.append(mut_aux)
            else:
                continue

    mut_def = [item for sublist in muts for item in ([sublist] if isinstance(sublist, list) else [sublist])]
    mut_df = pd.DataFrame(mut_def,columns=[
                "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                "Nt_mutation", "Mutation_frequency"
            ])
    return mut_df
