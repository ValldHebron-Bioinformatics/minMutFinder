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

import sys
import pandas as pd

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--tsv':
        in_tsv = sys.argv[i + 1]

SAMPLE = out_dir.split('/')[len(out_dir.split('/'))-1]
VARIANT_CALLING = out_dir + '/variant_calling'

ivar_tsv = pd.read_csv(in_tsv, sep='\t')

with open(VARIANT_CALLING + '/' + SAMPLE + '_prot_variants.vcf', "w") as vcf_file:
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
    vcf_file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for index, row in ivar_tsv.iterrows():
        chrom = row['REGION']
        pos = row['POS']
        ref = row['REF']
        alt = row['ALT']
        total_dp = row['TOTAL_DP']
        alt_freq = row['ALT_FREQ']

        info = f"DP={total_dp};AF={alt_freq}"

        vcf_file.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
