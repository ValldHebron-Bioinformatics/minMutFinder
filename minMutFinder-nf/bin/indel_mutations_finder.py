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


def indel_mutations_detector(vcf, ref, sample_id, protein_initial):
    mut = []
    ref_seq = ref
    df = vcf
    for index, row in df.iterrows():
        ref, alt, pos, allele_freq, allele_dp = df.REF.iloc[index], df.ALT.iloc[index], df.POS.iloc[index], df.AF.iloc[index], df.DP.iloc[index]

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
        aa_pos, protein, aux = int(pos // 3) + 1, protein_initial, 0

        if ((len(ref) > 1) or (len(alt) > 1) ) and (len(ref) != len(alt)):
            aux, POS = 1, pos
            if ((int((len(ref) - 1) % 3) == 0) and (int(len(ref) != 1))) or ((int((len(alt) - 1) % 3) == 0) and (int(len(alt) != 1))):

                if (int(len(ref)) > 3) or (int(len(alt)) > 3):
                    n_positional = n-1
                    if (len(ref) > len(alt)):
                        # ref_complete=ref_seq[POS-n_positional:POS+len(ref)+((len(ref)+1)%3)]
                        # alt_complete=ref_complete.replace(ref[1:],'')[0:3]
                        del_nt=ref[1:]
                        if n == 0:
                            ref_complete=del_nt
                            #alt_complete=ref_seq[POS+1:POS+len(ref)-1+4]
                        elif n == 1:
                            ref_complete=ref_seq[POS]+del_nt+ref_seq[POS+len(ref):POS+len(ref)+2]
                            #alt_complete=ref_seq[POS:POS+len(ref)-1+3]
                            #alt_complete=ref_seq[POS-1]+ins_nt+ref_seq[POS-1:POS+2]
                        else: # n == 2
                            ref_complete=ref_seq[POS-1:POS+1]+del_nt+ref_seq[POS+len(ref)]
                            #alt_complete=ref_seq[POS-1:POS+len(ref)-1+2]
                        alt_complete=ref_complete.replace(ref[1:],'')[0:3]
                        for i in range(0,len(ref_complete),3):
                            if len(ref_complete[i:i+3]) == 3:
                                prot_ref=translateDNA(ref_complete[i:i+3])
                                prot_sub=translateDNA(alt_complete)
                                if ((n != 0) and (i == 0)): # or (len(ref) > 4))  and (len(alt_complete) > 2):
                                    if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                                        aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                                    elif (aa_classes[prot_ref] == "Special case") and (aa_classes[prot_sub] == "Special case"):
                                        if (prot_ref != prot_sub):
                                            aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                                        else:
                                            aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

                                    else:
                                        aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

                                    if prot_ref == prot_sub:
                                        if ref_complete[i:i+3] != alt_complete:
                                            mut_aux = [
                                                sample_id, protein, "SYNONYMOUS",
                                                translateDNA(ref_complete[i:i+3]) + str(aa_pos) + translateDNA(alt_complete),
                                                aa_change, ref_complete[i:i+3] + str(pos+1) + alt_complete, allele_freq, allele_dp
                                            ]
                                            mut.append(mut_aux)
                                    else:
                                        mut_aux = [
                                            sample_id, protein, "NON_SYNONYMOUS",
                                            translateDNA(ref_complete[i:i+3]) + str(aa_pos) + translateDNA(alt_complete),
                                            aa_change, ref_complete[i:i+3] + str(pos+1) + alt_complete, allele_freq, allele_dp
                                        ]
                                        mut.append(mut_aux)                                        
                                else:
                                    if prot_ref != prot_sub or n == 0:
                                        
                                        mut_aux = [
                                            sample_id, protein, "NON_SYNONYMOUS",
                                            translateDNA(ref_complete[i:i+3]) + str((((POS-n_positional)+i)//3)+1) + '-',
                                            "Deletion", ref_complete[i:i+3] + str((POS-n_positional+1)+i) + '-', allele_freq, allele_dp
                                        ]
                                        mut.append(mut_aux)
                                aa_pos+=1

                    else:
                        ins_nt=alt[1:]
                        if n == 0:
                            ref_complete=ref_seq[POS+1:POS+4]
                            alt_complete=ins_nt
                        elif n == 1:
                            ref_complete=ref_seq[POS:POS+3]
                            alt_complete=ref_seq[POS]+ins_nt+ref_seq[POS+1:POS+3]
                            #alt_complete=ref_seq[POS-1]+ins_nt+ref_seq[POS-1:POS+2]
                        else: # n == 2
                            ref_complete=ref_seq[POS-1:POS+2]
                            alt_complete=ref_seq[POS-1:POS+1]+ins_nt+ref_seq[POS+1]
                            # alt_complete=ref_seq[POS-2]+ins_nt+ref_seq[POS-2:POS+1]
                        #ref_complete=ref_seq[POS-n_positional:POS+n+((n+1)%3)]

                        #alt_complete=ref_seq[POS-n_positional]+ins_nt+ref_seq[POS-n_positional+1:POS+n+((n+1)%3)]
                        #alt_complete=ref_seq[POS-n:POS+len(ref)+((len(ref)+1)%3)]
                        for i in range(0,len(alt_complete),3):
                            if len(alt_complete[i:i+3]) == 3:
                                prot_sub=translateDNA(alt_complete[i:i+3])
                                prot_ref=translateDNA(ref_complete)
                                if ((n != 0) and (i == 0)): # and (len(ref_complete) > 2):
                                    if (aa_classes[prot_ref] != aa_classes[prot_sub]):
                                        aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                                    elif (aa_classes[prot_ref] == "Special case") and (aa_classes[prot_sub] == "Special case"):
                                        if (prot_ref != prot_sub):
                                            aa_change = "Amino acid changed from " + aa_classes[prot_ref] + " to " + aa_classes[prot_sub]

                                        else:
                                            aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

                                    else:
                                        aa_change = "Amino acid did not change, it stayed " + aa_classes[prot_ref]

                                    if prot_ref == prot_sub:
                                        if alt_complete[i:i+3] != ref_complete:
                                            mut_aux = [
                                                sample_id, protein, "SYNONYMOUS",
                                                translateDNA(ref_complete) + str(aa_pos) + translateDNA(alt_complete[i:i+3]),
                                                aa_change, ref_complete + str(pos+1) + alt_complete[i:i+3], allele_freq, allele_dp
                                            ]
                                            mut.append(mut_aux)
                                    else:
                                        mut_aux = [
                                            sample_id, protein, "NON_SYNONYMOUS",
                                            translateDNA(ref_complete) + str(aa_pos) + translateDNA(alt_complete[i:i+3]),
                                            aa_change, ref_complete + str(pos+1) + alt_complete[i:i+3], allele_freq, allele_dp
                                        ]
                                        mut.append(mut_aux)                                        
                                else:
                                    if (prot_ref != prot_sub) or n == 0:
                                        mut_aux = [
                                            sample_id, protein, "NON_SYNONYMOUS",
                                            '-' + str((((POS-n_positional)+i)//3)+1) + translateDNA(alt_complete[i:i+3]),
                                            "Insertion", '-' + str((POS-n_positional+1)+i) + alt_complete[i:i+3], allele_freq, allele_dp
                                        ]
                                        mut.append(mut_aux)
                                aa_pos+=1

                else:
                    if (len(ref) > len(alt)):
                        mut_aux = [
                            sample_id, protein, "NON_SYNONYMOUS",
                            "del" + str(aa_pos), "Deletion",
                            ref + str(POS+1) + alt, allele_freq, allele_dp
                        ]
                        mut.append(mut_aux)

                    else:
                        mut_aux = [
                            sample_id, protein, "NON_SYNONYMOUS",
                            "ins" + str(aa_pos), "Insertion",
                            ref + str(POS+1) + alt, allele_freq, allele_dp
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
                    sample_id, protein, "NON_SYNONYMOUS",
                    "Frameshift" + str(aa_pos),
                    "Frameshift", ref + str(POS) + alt, allele_freq, allele_dp
                ]
                mut.append(mut_aux)
    mut_def = [item for sublist in mut for item in ([sublist] if isinstance(sublist, list) else [sublist])]
    mut_df = pd.DataFrame(mut_def,columns=[
                    "SampleID", "Protein", "Mutation_type", "Aa_change", "Amino_Acid_Property_Change ",
                    "Nt_mutation", "Mutation_frequency", "Mutation_depth"
                ])
    return mut_df
