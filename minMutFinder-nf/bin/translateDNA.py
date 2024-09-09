#!/usr/bin/env python3
import re

nt_degenerated = {
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT"
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


def translateDNA(seq):
    protein_sequence = ""
    # Generate protein sequence
    for i in range(0, len(seq)):
        if len(seq) < 3:
            protein_sequence += "-"
        elif re.search("N", seq[i:i + 3]):
            protein_sequence += "X"
        elif seq[i:i+3] == "---":
            protein_sequence += "-"
        elif re.search("-", seq[i:i+3]):
            protein_sequence += "?"
        elif re.search('[RYSWKMBDHV]', seq[i:i+3]):
            place = re.search('[RYSWKMBDHV]', seq[i:i+3]).start()
            aa_opts = []
            for v in nt_degenerated[seq[i+place:i+place+1]]:
                seq_degenerated = seq[i:i+3]
                seq_degenerated = seq_degenerated.replace(seq[i + place:i + place + 1], v)
                aa_opts.append(codon_code[seq_degenerated])
            for o in range(len(aa_opts)):
                if o == 0:
                    continue
                else:
                    # No codifican mismo aa
                    if aa_opts[o] != aa_opts[o-1]:
                        protein_sequence += "?"
                    # Si llega al final de la iteracion y no ha entrado en el bucle anterior, son todos iguales
                    else:
                        if o == len(aa_opts)-1:
                            protein_sequence += codon_code[seq_degenerated]
        else:
            if codon_code[seq[i:i+3]] == '*':
                return protein_sequence
            else:
                protein_sequence += codon_code[seq[i:i+3]]
        return protein_sequence
