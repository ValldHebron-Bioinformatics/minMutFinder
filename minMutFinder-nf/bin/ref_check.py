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
