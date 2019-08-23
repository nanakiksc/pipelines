#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import collections

chroms = collections.defaultdict(list)

with open(sys.argv[1]) as fin:
    for line in fin:
        chroms[line.split()[0]].append(line)

for chrom in sorted(chroms):
    for line in chroms[chrom]:
        sys.stdout.write(line)
