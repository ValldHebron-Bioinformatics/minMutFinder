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

import gzip
import sys
import os
import shutil
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from read_count import read_count

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--r1':
        r1 = sys.argv[i + 1]
    elif sys.argv[i] == '--r2':
        r2 = sys.argv[i + 1]
    elif sys.argv[i] == '--ref-seq':
        ref_seq = sys.argv[i + 1]
    elif sys.argv[i] == '--sample':
        SAMPLE = sys.argv[i + 1]
    elif sys.argv[i] == '--r1qc':
        FQ_R1 = sys.argv[i + 1]
    elif sys.argv[i] == '--r2qc':
        FQ_R2 = sys.argv[i + 1]
    elif sys.argv[i] == '--r1UPqc':
        FQ_R1_UP = sys.argv[i + 1]
    elif sys.argv[i] == '--r2UPqc':
        FQ_R2_UP = sys.argv[i + 1]


FASTQ = out_dir + '/fastq'
QC_DIR = out_dir + '/qc'
qc_metrics_file = QC_DIR + '/qc_metrics.csv'
PLOTS = out_dir + '/plots'
ASSEMBLY = out_dir + '/assembly'



R1_file, R2_file = SeqIO.parse(gzip.open(r1, 'rt'), 'fastq'), SeqIO.parse(gzip.open(r2, 'rt'), 'fastq')
R1_qc_paired_file = SeqIO.parse(gzip.open(FQ_R1, 'rt'), 'fastq')
R2_qc_paired_file = SeqIO.parse(gzip.open(FQ_R2, 'rt'), 'fastq')
R1_qc_file = SeqIO.parse(gzip.open(FQ_R1_UP, 'rt'), 'fastq')
R2_qc_file = SeqIO.parse(gzip.open(FQ_R2_UP, 'rt'), 'fastq')

reads_total = read_count(R1_file, R2_file)
reads_filt_qc_unpaired = read_count(R1_qc_file, R2_qc_file)
reads_filt_qc_paired = read_count(R1_qc_paired_file, R2_qc_paired_file)
reads_filt_qc = int(reads_filt_qc_paired) + int(reads_filt_qc_unpaired)

if os.path.isfile(QC_DIR + '/qc_metrics.csv') is False:
    with open(QC_DIR + '/qc_metrics.csv', 'w') as qc_metrics:
        print('QC metrics file created')
        qc_metrics.write('sample;test;score' + '\n')
        qc_metrics.write(SAMPLE + ';reads_total;' + str(reads_total) + '\n')
        qc_metrics.write(SAMPLE + ';reads_filt_qc;' + str(reads_filt_qc) + '\n')
        qc_metrics.write(SAMPLE + ';reads_filt_qc_paired;' + str(reads_filt_qc_paired) + '\n')
qc_metrics.close()

qc = pd.read_csv(qc_metrics_file, sep=';')

sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context("paper")
pl = sns.barplot(x=qc.test, y=qc.score)
pl.tick_params(axis='x', labelsize=5)
plt.savefig(PLOTS + '/qc_metrics.svg', format='svg', dpi=2400)
