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
