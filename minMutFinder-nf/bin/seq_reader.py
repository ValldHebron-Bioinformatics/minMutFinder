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


# Store sequence
def seq_reader(file):
    sequences, sequence_array, header = {}, [], ""
    with open(file, 'r') as infile:
        while True:
            line = infile.readline()
            if not line:
                seq = "".join(sequence_array)
                sequences[header] = seq
                break
            else:
                if line.startswith(">"):
                    if header == "":
                        header = line[1:].rstrip()
                    else:
                        seq = "".join(sequence_array)
                        sequences[header] = seq
                        header, sequence_array = line[1:].rstrip(), []
                else:
                    sequence_array.append(line.rstrip())
    infile.close()
    return sequences
