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

import pandas as pd

def expand_vcf_row(row):
    alts = row['ALT'].split(',')
    
    info_fields = {key: value for key, value in [field.split('=') for field in row['INFO'].split(';') if '=' in field]}
    
    afs = info_fields['AF'].split(',')
    dps = info_fields['DP'].split(',')
    
    if len(afs) == 1:
        afs = afs * len(alts)
    if len(dps) == 1:
        dps = dps * len(alts)

    expanded_rows = []
    for alt, af, dp in zip(alts, afs, dps):
        new_row = row.copy()
        new_row['ALT'] = alt

        info_fields['AF'] = af
        info_fields['DP'] = dp

        new_info = ';'.join([f"{key}={value}" for key, value in info_fields.items()])
        new_row['INFO'] = new_info
        
        expanded_rows.append(new_row)
    
    return expanded_rows
