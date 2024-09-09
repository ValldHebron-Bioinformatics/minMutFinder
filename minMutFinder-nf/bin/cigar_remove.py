#!/usr/bin/env python3


def cigar_remove(cigar_input, startH, new_cig, idx, number_class, H_exists, p):
    if 'H' in cigar_input:
        H_exists = True
        for c in p.finditer(cigar_input):
            if (number_class == 0) & (c.group() == "H"):
                print(cigar_input + '\nRemoving H flag...')
                cigar_input = cigar_input.replace('H', '', 1)
                endH = c.start()
                if ((endH - startH) > 2):
                    for i in range(0, len(cigar_input)):
                        if i not in range(startH + 1, endH):
                            new_cig += cigar_input[i]
                elif startH == 0:
                    for i in range(0, len(cigar_input)):
                        if i not in range(startH, endH):
                            new_cig += cigar_input[i]
                else:
                    for i in range(0, len(cigar_input)):
                        if i != idx:
                            new_cig += cigar_input[i]
                print(new_cig)
                return cigar_remove(new_cig, startH, '', 0, 0, H_exists, p)

                break
            else:
                startH = c.start()
            idx = c.start()+1
    else:
        if (H_exists):
            print('ALL Hs removed...')
            print(cigar_input)
        return cigar_input
