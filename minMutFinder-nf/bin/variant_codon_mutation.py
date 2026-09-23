#!/usr/bin/env python3

# This file is part of minMutFinder.
#
# minMutFinder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# minMutFinder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2024 Respiratory Viruses Unit, Microbiology Department, Vall d’Hebron Hospital Universitari, Vall d’Hebron Institut de Recerca (VHIR), Vall d’Hebron Barcelona Hospital Campus, Passeig Vall d’Hebron 119-129, 08035 Barcelona, Spain

import pandas as pd
from replacer_vcm import replacer_vcm
from all_muts_finder import all_muts_finder

aa_classes = {
    'A': 'Hydrophobic', 'V': 'Hydrophobic',
    'L': 'Hydrophobic', 'M': 'Hydrophobic',
    'I': 'Hydrophobic', 'S': 'Polar non charged',
    'T': 'Polar non charged', 'N': 'Polar non charged',
    'Q': 'Polar non charged', 'G': 'Special case',
    'C': 'Special case', 'P': 'Special case',
    'U': 'Special case', 'F': 'Aromatic Hydrophobic',
    'Y': 'Aromatic Hydrophobic', 'W': 'Aromatic Hydrophobic',
    'K': 'Positively charged', 'R': 'Positively charged',
    'H': 'Positively charged', 'D': 'Positively charged',
    'E': 'Positively charged', '*': 'Stop'
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


def variant_codon_mutation(codon_data, ref_seq, reads, aa_pos, sample_id, protein, AF, DP):
    muts = []

    # Data loading
    codon, fasta, aa_pos, FREQ_CUTOFF, DEPTH_CUTOFF = codon_data, ref_seq, int(aa_pos), float(AF), int(DP)

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
                percentage, protein_pos, dp = (n_reads_mut / count), protein, n_reads_mut
                if percentage >= FREQ_CUTOFF and dp >= DEPTH_CUTOFF:

                    protein_pos, mut_dic[mut], prot_sub = protein, percentage, codon_code[mut]
                    prot_ref, nt_mutation = codon_code[reference_codon], ''

                    mut_aa_dic[prot_sub] = prot_ref
                    for i in r:
                        if reference_codon[i] != mut[i]:
                            if nt_mutation == '':
                                nt_mutation += str(str(reference_codon[i]) + str(i + (aa_pos * 3) - 2) + str(mut[i]))
                            else:
                                nt_mutation += str(", " + str(reference_codon[i]) + str(i + (aa_pos * 3) - 2) + str(mut[i]))

                    if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                        aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                    elif (aa_classes[prot_ref] == "Special case") and (aa_classes[prot_sub] == "Special case"):
                        if (prot_ref != prot_sub):
                            aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                        else:
                            aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

                    else:
                        aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

                    if prot_ref != prot_sub:
                        mut_aux = [
                            sample_id, protein_pos, "NON_SYNONYMOUS",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage),
                            str(dp)
                        ]

                    else:
                        mut_aux = [
                            sample_id, protein_pos, "SYNONYMOUS",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage),
                            str(dp)
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
                fixed_mut_dp = float(codon.DP.iloc[index])
            else:
                r_not_fixed.append(AA_INDEX)
        protein_pos = protein

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
                            dp_mut = b["DP"].to_string(index=False).split('\n')[0]
                            dp_non_fixed = float(dp_mut)
                            fixed_mut_dp = fixed_mut_dp - dp_non_fixed
                            prot_sub = codon_code[mutated_codon]

                            if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                                aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                            elif (aa_classes[prot_ref] == "Special case") and (aa_classes[prot_sub] == "Special case"):
                                if (prot_ref != prot_sub):
                                    aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]
                                else:
                                    aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]
                            else:
                                aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

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
                            dp = dp_non_fixed
                            if percentage >= FREQ_CUTOFF and dp >= DEPTH_CUTOFF:
                                if prot_ref != prot_sub:
                                    mut_aux = [
                                        sample_id, protein_pos, "NON_SYNONYMOUS",
                                        prot_ref + str(aa_pos) + prot_sub,
                                        aa_change, nt_mutation, str(percentage),
                                        str(dp)
                                    ]
                                else:
                                    mut_aux = [
                                        sample_id, protein_pos, "SYNONYMOUS",
                                        prot_ref + str(aa_pos) + prot_sub,
                                        aa_change, nt_mutation, str(percentage),
                                        str(dp)
                                    ]
                                muts.append(mut_aux)
                else:
                    continue

        nt_fixed_list = []
        for i in r_fixed:
            prot_sub = codon_code[fixed_nts]
            if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]
            elif (aa_classes[prot_ref] == "Special case") and (aa_classes[prot_sub] == "Special case"):
                if (prot_ref != prot_sub):
                    aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]
                else:
                    aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]
            else:
                aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

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
                dp = fixed_mut_dp

                if percentage >= FREQ_CUTOFF and dp >= DEPTH_CUTOFF:
                    if prot_ref != prot_sub:
                        mut_aux = [
                            sample_id, protein_pos, "NON_SYNONYMOUS",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage),
                            str(dp)
                        ]
                    else:
                        mut_aux = [
                            sample_id, protein_pos, "SYNONYMOUS",
                            prot_ref + str(aa_pos) + prot_sub,
                            aa_change, nt_mutation, str(percentage),
                            str(dp)
                        ]
                    muts.append(mut_aux)
            else:
                continue

    mut_def = [item for sublist in muts for item in ([sublist] if isinstance(sublist, list) else [sublist])]
    mut_df = pd.DataFrame(mut_def,columns=[
                "SampleID", "Protein", "Mutation_type", "Aa_change", "Amino_Acid_Property_Change",
                "Nt_mutation", "Mutation_frequency", "Mutation_depth"
            ])
    return mut_df
