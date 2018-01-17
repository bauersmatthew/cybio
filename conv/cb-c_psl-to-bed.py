#!/usr/bin/env python

# Script: "cb-c_psl-to-bed.py"
# Requires: Python 3
# Purpose: To convert PSL files (e.g. from BLAT) to BED files. Bedops includes a
#          psl2bed utility (a.k.a. convert2bed), but it thinks that BED files
#          have 20 columns...
# Author: Matthew Bauer

import sys
from math import floor
import argparse

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
    ret['score'] = floor(1000*int(fields[0])/(
        int(fields[0])+int(fields[1]))) # "permille" matching

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

    # check the cutoff
    pm_match = i['score']
    pm_cutoff = floor(10*args.cutoff)
    if pm_match < pm_cutoff:
        return

    sys.stdout.write(
        ('{c}\t{s}\t{e}\t{n}\t{score}\t{r}'
         '\t{s}\t{e}\t0,0,0\t{bc}\t{bz}\t{bs}\n').format(
             c=i['chrom'], s=i['start'], e=i['end'], r=i['strand'],
             n=i['name'], bc=i['b_num'], bz=i['b_sizes'], bs=i['b_starts'],
             score=i['score']))

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

def valid_percent(f):
    """Raise an exception if f is not a valid percent."""
    f = float(f)
    if f < 0.0 or f > 100.0:
        raise RuntimeError('{} is not between 0 and 100!'.format(f))
    return f

arg_parser = argparse.ArgumentParser(
    description='Convert a PSL file to a BED file.',
    epilog=('Results are written to STDOUT. The score field is set to '
            'the \'permille\' match score of the alignment.'))
arg_parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.1.0')
arg_parser.add_argument('-i', '--input',
                        type=valid_file, metavar='in.psl',
                        help=('The input PSL file. If not specified, input is '
                              'read from STDIN.'))
arg_parser.add_argument('-c', '--cutoff',
                        type=valid_percent, metavar='%_matches', default=50.0,
                        help=('The cutoff for % matches in the alignment; '
                              'reads below this cutoff are discarded. '
                              'Default: 50.0'))
args = arg_parser.parse_args()

fin = sys.stdin
if args.input is not None:
    fin = open(args.input)

for line in fin:
    write_bed(parse_psl(line.strip()))
