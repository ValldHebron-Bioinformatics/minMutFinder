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
import re
from cigar_remove import cigar_remove

def reads_assembler(samfile, all_fastas):
    print('Beginning reads assembly of ' + str(samfile))
    infile = open(samfile, 'r')
    while infile:
        line = infile.readline()
        df_cigar = pd.DataFrame(columns=['NODE', 'INI', 'END', 'POS', 'CLASS', 'SEQ'])
        if (line == ""):
            break
        # Create segments df
        node, flag = line.split('\t')[0], line.split('\t')[1]
        if (int(flag) == -1):
            isit = 0
        else:
            isit = 1
            pos, pos_end = int(line.split('\t')[3]), int(line.split('\t')[3])
            cigar, sequence, p, idx, seq_init = line.split('\t')[5], line.split('\t')[9], re.compile("[a-zA-Z]"), 0, 0
            class_del, number_class, new_cig, startH, H_exists = 0, 0, '', 0, False

            if cigar == "*":
                continue

            cigar = cigar_remove(cigar, startH, new_cig, idx, number_class, H_exists, p)

            idx = 0
            for c in p.finditer(cigar):
                end, seq_end = int(pos) + int(cigar[idx:c.start()]) - 1, seq_init + int(cigar[idx:c.start()])
                if (number_class == 0) & (c.group() == "S"):
                    pos = int(pos) - int(cigar[idx:c.start()])

                if c.group() == "D":
                    seq, class_del, count_del = "X" * int(cigar[idx:c.start()]), 1, int(cigar[idx:c.start()])
                else:
                    if class_del == 1:
                        seq_init, seq_end = (seq_init - count_del), (seq_end - count_del)
                    seq, class_del = sequence[seq_init:seq_end], 0
                    # class_del = 0
                new_row = [node, pos, end, cigar[idx:c.start()], c.group(), seq]
                pos, idx, seq_init = (int(pos) + int(cigar[idx:c.start()])), (c.start() + 1), seq_end

                df_cigar.loc[len(df_cigar)] = new_row

            # 2. Process the data

            # print("Filtering data...")

            # Remove S and H categories (clipped sequences)
            df_segments = df_cigar.loc[(df_cigar['CLASS'] != 'S') & (df_cigar['CLASS'] != 'H'), ]
            # print(df_segments.sort_values('INI'))
            df_segments = df_segments.sort_values('INI')
            startN = df_segments['INI'].iloc[0]  # Correct the start of the consensus sequence to the first position defined
            df_segments['INI'], df_segments['END'] = (df_segments['INI'] - startN), (df_segments['END'] - startN)

            for f in range(0, len(df_segments)):
                if df_segments['CLASS'].iloc[f] == "I":
                    for r in range(f + 1, len(df_segments)):
                        if df_segments['NODE'].iloc[r] == df_segments['NODE'].iloc[f]:
                            continue
                        else:
                            df_segments['INI'].iloc[r] = df_segments['INI'].iloc[r] + int(df_segments['POS'].iloc[f])
                            df_segments['END'].iloc[r] = df_segments['END'].iloc[r] + int(df_segments['POS'].iloc[f])

            # 3. Generate consensus sequence
            # print("Generating consensus sequence...")
            df_consensus = pd.DataFrame(columns=['POS', 'NT'])

            for i in range(0, len(df_segments)):
                for j in range(df_segments['INI'].iloc[i], df_segments['END'].iloc[i] + 1):
                    if 0 <= j - df_segments['INI'].iloc[i] < len(df_segments['SEQ'].iloc[i]):
                        new_row_c = [j, df_segments['SEQ'].iloc[i][j - df_segments['INI'].iloc[i]]]
                        df_consensus.loc[len(df_consensus)] = new_row_c

            # Add 'N' in the positions that are not defined.
            for p in range(df_consensus['POS'].iloc[0], df_consensus['POS'].iloc[-1]):
                if p not in df_consensus['POS'].tolist():
                    new_row_c = [p, 'N']
                    df_consensus.loc[len(df_consensus)] = new_row_c

            df_consensus = df_consensus.drop_duplicates(subset=['POS'])
            df_consensus = df_consensus.sort_values('POS')
            consensus_list = df_consensus['NT'].tolist()
            consensus_seq = ''.join(consensus_list)

        # with open(out_file, 'w') as outfile:
        spaces = "-" * (pos_end - 1)
        ending = "-" * (3 - (len(spaces + consensus_seq) % 3))
        fasta_seq = '>' + node + '\n' + spaces + consensus_seq + ending + '\n'
        if (isit == 1):
            all_fastas += fasta_seq
    return all_fastas
