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
from translateDNA import translateDNA

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


def indel_mutations_detector(vcf, ref, sample_id, gene_initial):
    mut = []
    ref_seq = ref
    df = vcf
    for index, row in df.iterrows():
        ref, alt, pos, allele_freq = df.REF.iloc[index], df.ALT.iloc[index], df.POS.iloc[index], df.AF.iloc[index]
        #alt, ref, pos, allele_freq = df.REF.iloc[index], df.ALT.iloc[index], df.POS.iloc[index], df.AF.iloc[index]

        n = int(pos) % 3

        if len(ref) > 1:
            if n == 1:
                ALT = str(alt) + str(ref_seq[int(pos) + len(ref) - 1:int(pos) + len(ref) + 1])

            elif n == 2:
                ALT = str(ref_seq[int(pos) - 1]) + str(alt) + str(ref_seq[int(pos) + len(ref) - 1])

            else:
                ALT = str(ref_seq[int(pos) - 2:int(pos)]) + str(alt)

        else:
            ALT = alt

        pos = int(pos) - 1
        aa_pos, gene, aux = int(pos // 3) + 1, gene_initial, 0

        if ((len(ref) > 1) or (len(alt) > 1) ) and (len(ref) != len(alt)):
            aux, POS = 1, pos
            if ((int((len(ref) - 1) % 3) == 0) and (int(len(ref) != 1))) or ((int((len(alt) - 1) % 3) == 0) and (int(len(alt) != 1))):

                if (int(len(ref)) > 3) or (int(len(alt)) > 3):
                    n_positional = n-1
                    if (len(ref) > len(alt)):
                        ref_complete=ref_seq[POS-n_positional:POS+len(ref)+((len(ref)+1)%3)]
                        alt_complete=ref_complete.replace(ref[1:],'')[0:3]
                        for i in range(0,len(ref_complete),3):
                            if len(ref_complete[i:i+3]) == 3:
                                prot_ref=translateDNA(ref_complete[i:i+3])
                                prot_sub=translateDNA(alt_complete)
                                if (((n != 0) or (len(ref) > 4)) and (i == 0)) and (len(alt_complete) > 2):
                                    if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                                        aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                                    elif (aa_classes[prot_ref] == "Special_case") and (aa_classes[prot_sub] == "Special_case"):
                                        if (prot_ref != prot_sub):
                                            aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                                        else:
                                            aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                                    else:
                                        aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                                    if prot_ref == prot_sub:
                                        if ref_complete[i:i+3] != alt_complete:
                                            mut_aux = [
                                                sample_id, gene, "NO_MUTATION",
                                                translateDNA(ref_complete[i:i+3]) + str(aa_pos) + translateDNA(alt_complete),
                                                aa_change, ref_complete[i:i+3] + str(pos+1) + alt_complete, allele_freq
                                            ]
                                            mut.append(mut_aux)
                                    else:
                                        mut_aux = [
                                            sample_id, gene, "MUTATION",
                                            translateDNA(ref_complete[i:i+3]) + str(aa_pos) + translateDNA(alt_complete),
                                            aa_change, ref_complete[i:i+3] + str(pos+1) + alt_complete, allele_freq
                                        ]
                                        mut.append(mut_aux)                                        
                                else:
                                    if prot_ref != prot_sub:
                                        
                                        mut_aux = [
                                            sample_id, gene, "MUTATION",
                                            translateDNA(ref_complete[i:i+3]) + str((((POS-n_positional)+i)//3)+1) + '-',
                                            "Deletion", ref_complete[i:i+3] + str((POS-n_positional+1)+i) + '-', allele_freq
                                        ]
                                        mut.append(mut_aux)
                                aa_pos+=1

                    else:
                        alt_complete=ref_seq[POS-n:POS+len(ref)+((len(ref)+1)%3)]
                        ref_complete=alt_complete.replace(alt[1:],'')[0:3]
                        for i in range(0,len(alt_complete),3):
                            if len(alt_complete[i:i+3]) == 3:
                                prot_sub=translateDNA(alt_complete[i:i+3])
                                prot_ref=translateDNA(ref_complete)
                                if (((n != 0) or (len(alt) > 4)) and (i == 0)) and (len(ref_complete) > 2):
                                    if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                                        aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                                    elif (aa_classes[prot_ref] == "Special_case") and (aa_classes[prot_sub] == "Special_case"):
                                        if (prot_ref != prot_sub):
                                            aa_change = aa_classes[prot_ref] + " --> " + aa_classes[prot_sub]

                                        else:
                                            aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                                    else:
                                        aa_change = "No_Change, stayed " + aa_classes[prot_ref]

                                    if prot_ref == prot_sub:
                                        if alt_complete[i:i+3] != ref_complete:
                                            mut_aux = [
                                                sample_id, gene, "NO_MUTATION",
                                                translateDNA(ref_complete) + str(aa_pos) + translateDNA(alt_complete[i:i+3]),
                                                aa_change, ref_complete + str(pos+1) + alt_complete[i:i+3], allele_freq
                                            ]
                                            mut.append(mut_aux)
                                    else:
                                        mut_aux = [
                                            sample_id, gene, "MUTATION",
                                            translateDNA(ref_complete) + str(aa_pos) + translateDNA(alt_complete[i:i+3]),
                                            aa_change, ref_complete + str(pos+1) + alt_complete[i:i+3], allele_freq
                                        ]
                                        mut.append(mut_aux)                                        
                                else:
                                    if prot_ref != prot_sub:
                                        mut_aux = [
                                            sample_id, gene, "MUTATION",
                                            '-' + str((((POS-n_positional)+i)//3)+1) + translateDNA(alt_complete[i:i+3]),
                                            "Insertion", '-' + str((POS-n_positional+1)+i) + alt_complete[i:i+3], allele_freq
                                        ]
                                        mut.append(mut_aux)
                                aa_pos+=1

                else:
                    if (len(ref) > len(alt)):
                        mut_aux = [
                            sample_id, gene, "MUTATION",
                            "del" + str(aa_pos), "Deletion",
                            ref + str(POS+1) + alt, allele_freq
                        ]
                        mut.append(mut_aux)

                    else:
                        mut_aux = [
                            sample_id, gene, "MUTATION",
                            "ins" + str(aa_pos), "Insertion",
                            ref + str(POS+1) + alt, allele_freq
                        ]
                        mut.append(mut_aux)
            else:

                if (len(ref) > len(alt)):
                    longest=ref
                    if n == 2:
                        prot = str(ref_seq[pos - 1]) + str(ALT) + str(ref_seq[pos + 1])
                        ref_prot = str(ref_seq[pos - 1:pos + 2])

                        prot_sub = translateDNA(prot)
                        prot_ref, POS = translateDNA(ref_prot), int(pos) + 2

                    elif n == 0:
                        prot = str(ref_seq[pos - 2:pos]) + str(ALT)
                        ref_prot = str(ref_seq[pos - 2:pos + 1])

                        prot_sub = translateDNA(prot)
                        prot_ref, POS = translateDNA(ref_prot), int(pos) + 1

                    else:
                        prot = str(ALT) + str(ref_seq[pos + 1:pos + 3])
                        ref_prot = str(ref_seq[pos:pos + 3])

                        prot_sub = translateDNA(prot)
                        prot_ref, POS = translateDNA(ref_prot), int(pos) + 3
                else:
                    longest=alt
                    if n == 2:
                        ref_prot = str(ref_seq[pos - 1]) + str(ref) + str(ref_seq[pos + 1])
                        prot = str(ref_seq[pos - 1:pos + 2])
    
                        prot_sub = translateDNA(prot)
                        prot_ref, POS = translateDNA(ref_prot), int(pos) + 2
    
                    elif n == 0:
                        ref_prot = str(ref_seq[pos - 2:pos]) + str(ref)
                        prot = str(ref_seq[pos - 2:pos + 1])
    
                        prot_sub = translateDNA(prot)
                        prot_ref, POS = translateDNA(ref_prot), int(pos) + 1
    
                    else:
                        ref_prot = str(ref) + str(ref_seq[pos + 1:pos + 3])
                        prot = str(ref_seq[pos:pos + 3])
    
                        prot_sub = translateDNA(prot)
                        prot_ref, POS = translateDNA(ref_prot), int(pos) + 3

                mut_aux = [
                    sample_id, gene, "MUTATION",
                    "Frameshift" + str(aa_pos),
                    "Frameshift", ref + str(POS) + alt, allele_freq
                ]
                mut.append(mut_aux)

    mut_def = [item for sublist in mut for item in ([sublist] if isinstance(sublist, list) else [sublist])]
    mut_df = pd.DataFrame(mut_def,columns=[
                    "SampleID", "Gene", "Mutation_type", "Aa_change", "Type_of_aa_change ",
                    "Nt_mutation", "Mutation_frequency"
                ])
    return mut_df
