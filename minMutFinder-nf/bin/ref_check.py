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
import sys

for i in range(len(sys.argv)):
    if sys.argv[i] == '--ref-seq':
        ref_seq = sys.argv[i + 1]

print('Checking reference sequence file...')
if ('.fasta' not in ref_seq) or ('.fas' not in ref_seq):
    print('Reference protein sequence in incorrect format. \nIt has to be in .fasta or .fas...')
    sys.exit()
else:
    print('Reference protein ' + ref_seq + ' sequence in correct format.\nContinuing...')
