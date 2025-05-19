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
# Copyright (C) 2024 Ignasi Prats Méndez

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


def minority_mutations_detector(vcf, ref, sample_id, gene_initial):
    mut = []
    ref_seq = ref
    df = vcf
    for index, row in df.iterrows():

        ref, alt, pos, allele_freq = df.REF.iloc[index], df.ALT.iloc[index], df.POS.iloc[index], df.AF.iloc[index]

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


        if (n == 1) and (aux <= 0):
            n = 0  # first nucleotide

            prot = str(ALT) + str(ref_seq[pos + 1:pos + 3])

            prot_sub = codon_code[prot]

            ref_prot = str(ref_seq[pos:pos + 3])
            prot_ref = codon_code[ref_prot]
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
                mut = [
                    sample_id, gene, "MUTATION",
                    prot_ref + str(aa_pos) + prot_sub,
                    aa_change, ref + str(pos+1) + alt, allele_freq
                ]
            else:
                mut = [
                    sample_id, gene, "NO_MUTATION",
                    prot_ref + str(aa_pos) + prot_sub,
                    aa_change, ref + str(pos+1) + alt, allele_freq
                ]

        elif (n == 0) and (aux <= 0):
            n = 3  # last nucleotide

            prot = str(ref_seq[pos - 2:pos]) + str(ALT)

            prot_sub = codon_code[prot]

            ref_prot = str(ref_seq[pos - 2:pos + 1])
            prot_ref = codon_code[ref_prot]
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
                mut = [
                    sample_id, gene, "MUTATION",
                    prot_ref + str(aa_pos) + prot_sub,
                    aa_change, ref + str(pos+1) + alt, allele_freq
                ]
            else:
                mut = [
                    sample_id, gene, "NO_MUTATION",
                    prot_ref + str(aa_pos) + prot_sub,
                    aa_change, ref + str(pos+1) + alt, allele_freq
                ]

        else:
            if aux <= 0:
                n = 2  # nucleotide in the middle
                prot = str(ref_seq[pos -1 ]) + str(ALT) + str(ref_seq[pos + 1])
                prot_sub = codon_code[prot]

                ref_prot = str(ref_seq[pos - 1:pos + 2])
                prot_ref = codon_code[ref_prot]
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
                    mut = [
                        sample_id, gene, "MUTATION",
                        prot_ref + str(aa_pos) + prot_sub,
                        aa_change, ref + str(pos+1) + alt, allele_freq
                    ]
                else:
                    mut = [
                        sample_id, gene, "NO_MUTATION",
                        prot_ref + str(aa_pos) + prot_sub,
                        aa_change, ref + str(pos+1) + alt, allele_freq
                    ]
    return mut
