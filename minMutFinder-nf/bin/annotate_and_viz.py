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

import os
import pandas as pd
import sys
from visualizer import visualizer
from seq_reader import seq_reader

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--annotate':
        r_file = sys.argv[i + 1]
    elif sys.argv[i] == '--ref-seq':
        ref_seq = sys.argv[i + 1]
    elif sys.argv[i] == '--prot-names':
        names = sys.argv[i + 1]
    elif sys.argv[i] == '--syn-muts':
        syn_muts = sys.argv[i + 1]

SAMPLE = out_dir.split('/')[len(out_dir.split('/'))-1]
MUTATIONS = out_dir + '/mutations'
ASSEMBLY = out_dir + '/assembly'
QC_DIR = out_dir + '/qc'
qc_metrics_file = QC_DIR + '/qc_metrics.csv'
PLOTS = out_dir + '/plots'
PROT_REF = seq_reader(ref_seq)
NAMES = pd.read_csv(names,header=None)
names_list = list(NAMES[0])

lengths = []
c_row=1
for index in NAMES.index:
    if len(PROT_REF[NAMES.iloc[index].item()]) %3 == 0:
        end=len(PROT_REF[NAMES.iloc[index].item()])//3
    else:
        end=(len(PROT_REF[NAMES.iloc[index].item()])//3)+1

    lengths.append(list([NAMES.iloc[index].item(), c_row, 0, end]))
    c_row += 1

lengths_df = pd.DataFrame(lengths,columns=['ref','row','init','end'])

if os.path.isfile(ASSEMBLY + '/' + SAMPLE + '_depth_consensus.tsv'):
    df_depth_check = True
#     df_depth = pd.read_csv(ASSEMBLY + '/' + SAMPLE + '_depth_consensus.tsv', sep='\t')
else:
    df_depth_check = False
#     df_depth = False

df_muts_ALL = pd.read_csv(MUTATIONS + '/' + SAMPLE + '_mutations.csv', sep=';')

if syn_muts == 'no':
    df_muts_ALL = df_muts_ALL[df_muts_ALL.Mutation_type != "NO_MUTATION"]


annotate = pd.read_csv(r_file, sep='\t')
resistant_df = pd.DataFrame()
resistant_df['Resistant'] = ''
resistant_df['Resistant_mutation'] = ''
more_than_one = []
one = []
for index in annotate.index:
    row = annotate.iloc[index].str.split('\t')
    dfr = pd.DataFrame(row)
    dfr = dfr.T
    dfr.columns = annotate.columns
    muts = dfr['mutation']
    if len(str(muts).split('+')) > 1:
        more_than_one.append(list(annotate.iloc[index].T))
    else:
        one.append(list(annotate.iloc[index].T))

if (one):
    one = pd.DataFrame(one)
    one.columns = annotate.columns
    print(one)
    for index in one.index:
        print(df_muts_ALL.Gene)
        if (
            df_muts_ALL.Gene[df_muts_ALL.Gene == str(one.name[index])].shape[0] > 0
        ) and (
            df_muts_ALL[df_muts_ALL.Aa_change == str(one.mutation[index])].shape[0] > 0
        ):
            row = df_muts_ALL[df_muts_ALL.Aa_change == str(one.mutation[index])].reset_index(drop=True)
            row.columns = df_muts_ALL.columns
            row['Resistant'] = ['yes']
            row['Resistant_mutation'] = [str(one.mutation[index])]
            if resistant_df.size > 0:
                resistant_df = pd.concat([resistant_df, row])
            else:
                resistant_df = row

if (more_than_one):
    more_than_one = pd.DataFrame(more_than_one)
    more_than_one.columns = annotate.columns
    for index in more_than_one.index:
        all_muts = False
        for muts in more_than_one.mutation[index].split('+'):
            if (
                df_muts_ALL.Gene[df_muts_ALL.Gene == more_than_one.name[index]].shape[0] > 0
            ) and (
                df_muts_ALL[df_muts_ALL.Aa_change == muts].shape[0] > 0
            ):
                all_muts = True
            else:
                all_muts = False
                continue

        if (all_muts):
            for muts in more_than_one.mutation[index].split('+'):
                if (
                    df_muts_ALL.Gene[df_muts_ALL.Gene == more_than_one.name[index]].shape[0] > 0
                ) and (
                    df_muts_ALL[df_muts_ALL.Aa_change == muts].shape[0] > 0
                ):
                    row = df_muts_ALL[df_muts_ALL.Aa_change == muts].reset_index(drop=True)
                    row.columns = df_muts_ALL.columns
                    row['Resistant'] = ['yes']
                    row['Resistant_mutation'] = more_than_one.mutation[index]
                    if resistant_df.size > 0:
                        resistant_df = pd.concat([resistant_df, row])
                    else:
                        resistant_df = row

resistant_df = resistant_df.reset_index(drop=True)
print('Resistant mutations found:')
print(resistant_df)
# PLOTS
# visualizer(df_depth_check, df_depth, qc_metrics_file, df_muts_ALL, PLOTS, PROT_REF.keys(), resistant_df)
visualizer(df_depth_check, qc_metrics_file, df_muts_ALL, PLOTS, MUTATIONS, names_list, resistant_df, SAMPLE,syn_muts, lengths_df)

