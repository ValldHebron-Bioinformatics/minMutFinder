#!/usr/bin/env python3


def read_count(R1f, R2f):
    r1r = 0
    for rr1 in R1f:
        r1r += 1

    r2r = 0
    for rr2 in R2f:
        r2r += 1

    rt = r1r + r2r
    return int(rt)
