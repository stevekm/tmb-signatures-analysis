#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Count the total number of bases covered in a .bed file
"""
import sys
bedfile = sys.argv[1]
total = 0
with open(bedfile) as f:
    for line in f:
        parts = line.split()
        start = int(parts[1])
        stop = int(parts[2])
        dist = stop - start
        total = total + abs(dist)
print(total)
