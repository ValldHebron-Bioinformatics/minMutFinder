#!/usr/bin/env python3

from replacer_vcm import replacer_vcm


def all_muts_finder(alone_mut, r, codon_dic, muts_list):
    for nt in r:
        for c in codon_dic[nt]:
            mutated_codon = replacer_vcm(alone_mut, c, nt)
            if "Series" in mutated_codon:
                continue
            else:
                if mutated_codon not in muts_list:
                    muts_list.append(mutated_codon)
                alone_mut_aux = mutated_codon
                for nt1 in r:
                    for c1 in codon_dic[nt1]:
                        mutated_codon = replacer_vcm(alone_mut_aux, c1, nt1)
                        if "Series" in mutated_codon:
                            continue
                        else:
                            if mutated_codon not in muts_list:
                                muts_list.append(mutated_codon)
