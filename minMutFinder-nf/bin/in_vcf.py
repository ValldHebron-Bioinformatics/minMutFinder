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

import shutil
import sys

for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--vcf':
        in_vcf = sys.argv[i + 1]

SAMPLE = out_dir.split('/')[len(out_dir.split('/'))-1]
VARIANT_CALLING = out_dir + '/variant_calling'

if (in_vcf != VARIANT_CALLING + '/' + SAMPLE + '_prot_variants.vcf'):
    shutil.copyfile(in_vcf, VARIANT_CALLING + '/' + SAMPLE + '_prot_variants.vcf')
