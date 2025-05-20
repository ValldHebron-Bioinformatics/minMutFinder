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

import os
import sys
import shutil
import pysam
import pandas as pd

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--areads':
        areads = sys.argv[i + 1]

SAMPLE = out_dir.split('/')[len(out_dir.split('/'))-1]
QC_DIR = out_dir + '/qc'
qc_metrics_file = QC_DIR + '/QC_metrics.csv'
VARIANT_CALLING = out_dir + '/variant_calling'

if ".sam" in areads:
    print(areads + ' In sam format')
    if os.path.isfile(VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.sam'):
        print(VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.sam already exists, continuing...')
    else:
        shutil.copyfile(areads, VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.sam')
    pysam.view('-h', '-o', VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.bam', areads, catch_stdout=False)

elif ".bam" in areads: 
    print(areads + ' In bam format')
    if os.path.isfile(VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.bam'):
        print(VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.sam already exists, continuing...')
    else:
        shutil.copyfile(areads, VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.bam')
    pysam.view('-h', '-o', VARIANT_CALLING + '/' + SAMPLE + '_indelqual_prot.sam', areads, catch_stdout=False)

else:
    print("Incorrect file given!")
