#!/usr/bin/env python

# Script: "11-16_get-all-dscam-exons.py"
#         (specific ~one-time use)
# Requires: Python 3
# Purpose: Sensibly extract a BED record per Dscam1 exon in the 75 Dscam1
#          transcripts found in the Ensembl Dme v89 annotation.
# Usage: 11-16_get-all-dscam-exons.py Drosophila_melanogaster.BDGP6.89.gtf
# Author: Matthew Bauer

import sys

def process_gtf_record(rec):
    """Extract relevant information from the given GTF record."""
    fields = rec.split('\t')
    info = {}
    info['chr'] = fields[0]
    info['feature'] = fields[2]
    info['start'] = int(fields[3])-1 # conv to 0-based
    info['end'] = int(fields[4]) # no conversion needed
    for attr in fields[8].split(';'):
        parts = attr.split(' ')
        info[parts[0]] = parts[1].strip('"')
    return info

exons = [[]]*25 # DO NOT USE THE FIRST ELEMENT (exons[0])
# extract info from file
with open(sys.argv[1]) as fin:
    for line in fin:
        info = process_gtf_record(line)
        if (info['gene_name'] == 'Dscam1' and info['feature'] == 'exon'):
            interval = (info['start'], info['end'])
            exon_num = int(info['exon_number'])
            if interval not in exons[exon_num]:
                exons[exon_num].append(interval)
# sort all the intervals
for intvls in exons:
    intvls.sort(key=lambda i: i[0])
# write BED records for them
for exon_num, exon_vers in enumerate(exons):
    for ver_num, ver in exon_vers:
        sys.stdout.write(
            '2R\t{start}\t{end}\t{exon}.{ver}\n'.format(
                start=ver[0], end=ver[1], exon=exon_num, ver=ver_num))
