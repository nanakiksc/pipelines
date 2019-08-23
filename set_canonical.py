#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Create a subset of annovar's ensGene table
# containing only canonical transcripts
# as described in knownCanonical_hg38.bed

import sys
import gzip

rp = set(['RP1-', 'RP3-', 'RP4-', 'RP5-', 'RP6-', 'RP11-', 'RP13-'])

canonicals = set()
with gzip.open('knownCanonical_hg38.bed.gz') as kc:
    for i, line in enumerate(kc):
        canonicals.add(line.split()[4])
assert len(canonicals) == i + 1

with open(sys.argv[1]) as eg:
    for line in eg:
        sline = line.split()
        if sline[12][:4] in rp or sline[12][:5] in rp:
            continue
        if sline[1] in canonicals:
            sys.stdout.write(line)
