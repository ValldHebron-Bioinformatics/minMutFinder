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

#!/usr/bin/env python3
import os
import pandas as pd
import sys
from visualizer import visualizer
from seq_reader import seq_reader

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
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
NAMES = pd.read_csv(names,header=None)
names_list = list(NAMES[0])
PROT_REF = seq_reader(ref_seq)
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
    # df_depth = pd.read_csv(ASSEMBLY + '/' + SAMPLE + '_depth_consensus.tsv', sep='\t')
else:
    df_depth_check = False
    # df_depth = False

annotate_df = pd.DataFrame()
df_muts_ALL = pd.read_csv(MUTATIONS + '/' + SAMPLE + '_mutations.csv', sep=';')

if syn_muts == 'no':
    df_muts_ALL = df_muts_ALL[df_muts_ALL.Mutation_type != "NO_MUTATION"]

# visualizer(df_depth_check, df_depth, qc_metrics_file, df_muts_ALL, PLOTS, PROT_REF.keys(), annotate_df)
visualizer(df_depth_check, qc_metrics_file, df_muts_ALL, PLOTS, MUTATIONS, names_list, annotate_df, SAMPLE, syn_muts, lengths_df)
