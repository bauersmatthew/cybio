#!/usr/bin/env python

# Script: "cb-u-bed_set-strand.py"
# Requires: Python 3
# Purpose: This script sets the strand of each record in the given BED file
#          to be the same. This is useful for prepping liftOver, which really
#          messes up - strand stuff.
# Author: Matthew Bauer

import argparse
import os
import sys

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

def valid_strand(c):
    """Raise an exception if c is not a valid strand."""
    if c not in ('+', '-'):
        raise RuntimeError("'{}' is not a valid strand!".format(c))
    return c

arg_parser = argparse.ArgumentParser(
    description='Set the strand of every record to be the same.',
    epilog='Results are written to STDOUT.')
arg_parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0.0')
arg_parser.add_argument('-i', '--input', type=valid_file, required=False,
                        metavar='in.bed',
                        help=('The input BED file. If not given, STDIN '
                              'is used.'))
arg_parser.add_argument('strand', type=valid_strand,
                        help='+ or -')
args = arg_parser.parse_args()

fin = sys.stdin
if args.input is not None:
    fin = open(args.input)

for line in fin:
    f = line.rstrip().split()
    f[5] = args.strand
    sys.stdout.write('{}\n'.format('\t'.join(f)))

if args.input is not None:
    fin.close()
