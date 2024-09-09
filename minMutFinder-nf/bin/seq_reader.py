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
