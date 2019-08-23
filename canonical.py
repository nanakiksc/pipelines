#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Create a subset of annovar's refSeq table
# containing only canonical transcripts
# as described in NM.txt

import sys
from collections import defaultdict

nm_dict = {} # geneSymbol to refSeq dict.
with open('NM_full.txt') as nm:
    for line in nm:
        sline = line.rstrip().split(':')
        gene_symbol = sline[0]
        ref_seq = sline[1]
        nm_dict[gene_symbol] = ref_seq.split('.')[0] # Remove version number.

with open(sys.argv[1]) as rs:
    for line in rs:
        sline = line.split('\t')
        gene_symbol = sline[12]
        ref_seq = sline[1]

        if gene_symbol in nm_dict and ref_seq == nm_dict[gene_symbol]:
            print line.rstrip()
