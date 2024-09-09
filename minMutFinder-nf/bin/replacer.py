#!/usr/bin/env python3


def replacer(s, newstring, index):
    if s[index:index+1:] != "X":
        # if not erroring, but the index is still not in the correct range..
        if index < 0:  # add it to the beginning
            return newstring + s
        if index > len(s):  # add it to the end
            return s + newstring
        # insert the new string between "slices" of the original
        return (s[:index] + newstring + s[index + 1:])
