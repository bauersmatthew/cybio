#!/usr/bin/env python

# Script: "cb-c_psl-to-bed.py"
# Requires: Python 3
# Purpose: To convert PSL files (e.g. from BLAT) to BED files. Bedops includes a
#          psl2bed utility (a.k.a. convert2bed), but it thinks that BED files
#          have 20 columns...
# Author: Matthew Bauer

import sys
from math import floor

def parse_psl(line):
    """Parse one psl record. Return a dictionary of attributes."""
    if not line:
        return None
    ret = {}
    fields = line.split()
    if len(fields) != 21:
        raise RuntimeError('Invalid PSL line!')
    ret['name'] = fields[9].strip()
    ret['strand'] = fields[8].strip()[-1]
    ret['chrom'] = fields[13].strip()
    ret['start'] = fields[15].strip()
    ret['end'] = fields[16].strip()
    ret['b_num'] = fields[17].strip()
    ret['score'] = ','.join(fields[:2]) # matches,mismatches

    # block coordinates are weird; if on negative strand they are flipped
    l_sz = fields[18].strip().strip(',').split(',')
    if ret['strand'] == '-':
        l_sz = list(reversed(l_sz))
    ret['b_sizes'] = ','.join(l_sz)

    l_st = []
    s = int(ret['start'])
    for x in fields[20].strip().strip(',').split(','):
        x_fix = int(x)
        if ret['strand'] == '-':
            x_fix = int(fields[14].strip()) - int(x)
        x_fix -= s
        l_st.append(str(x_fix))
    if ret['strand'] == '-':
        l_st = list(reversed(l_st))
    if ret['strand'] == '-':
        for i, v in enumerate(l_st):
            l_st[i] = str(int(v)-int(l_sz[i]))
    ret['b_starts'] = ','.join(l_st)

    return ret

def write_bed(info):
    """Write a BED line to STDOUT given a dictionary of attributes."""
    if info is None:
        return
    i = info
    sys.stdout.write(
        ('{c}\t{s}\t{e}\t{n}\t{score}\t{r}'
         '\t{s}\t{e}\t0,0,0\t{bc}\t{bz}\t{bs}\n').format(
             c=i['chrom'], s=i['start'], e=i['end'], r=i['strand'],
             n=i['name'], bc=i['b_num'], bz=i['b_sizes'], bs=i['b_starts'],
             score=i['score']))

if len(sys.argv) > 1:
    sys.stderr.write('Usage: cb-c_psl-to-bed.py < in.psl > out.bed\n')
    sys.exit(0)

for line in sys.stdin:
    write_bed(parse_psl(line.strip()))
