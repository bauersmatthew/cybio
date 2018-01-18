#!/usr/bin/env python

# Script: "cb-np_filter-isoforms.py"
# Requires: Python 3
# Purpose: This script, given the output of cb-np_count-isoforms.py, will remove
#          unwanted isoforms.
# Author: Matthew Bauer

import argparse
import sys
import os.path
import re

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

def match_in(l, q):
    """Check if there exists a match in l to REGEX q."""
    for x in l:
        m = re.match(q, x)
        if m is not None and m.start() == 0 and m.end() == len(x):
            return True
    return False

arg_parser = argparse.ArgumentParser(
    description='Filter a table of isoforms.',
    epilog='Results are written to STDOUT.')
arg_parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0.0')
arg_parser.add_argument('isoforms',
                        type=valid_file,
                        help='The output file from cb-np_count-isoforms.py')
arg_parser.add_argument('filter',
                        type=str,
                        help=('The filter to be used. Format: \'A;B;C;...\' '
                              'where each component is \'[+-]feature\'. For '
                              'example: \'+11.4;-13.1\'.'))
arg_parser.add_argument('-n', '--negate',
                        action='store_true', default=False,
                        help=('Keep isoforms that do NOT match to to the '
                              'filter. Default: keep isoforms that do.'))
arg_parser.add_argument('-d', '--delim',
                        type=str, default=';',
                        help=('Set the delimiter that should be usd to '
                              'separate features in the filter string. '
                              'Default: ;'))
arg_parser.add_argument('-r', '--regex',
                        action='store_true', default=False,
                        help='Interpret the filter components as REGEXes.')
args = arg_parser.parse_args()

args.filter = args.filter.split(args.delim)
for f in args.filter:
    if len(f) < 2 or f[0] not in ('+', '-'):
        sys.stderr.write('Invalid filter string!\n')
        sys.exit(-1)

def filter_one(line):
    line = line.rstrip()
    if not line:
        return ''
    isos = line.split('\t')[0].split(',')
    for f in args.filter:
        s = f[0]
        n = f[1:]
        if args.negate:
            s = '+' if s == '-' else '-'
        found = None
        if not args.regex:
            found = n in isos
        else:
            found = match_in(isos, n)
        if (s == '+' and not found) or (s == '-' and found):
            # failed!
            return ''
    return line + '\n'

# read isoform table
iso_counts = []
with open(args.isoforms) as fin:
    for line in fin:
        sys.stdout.write(filter_one(line))

sys.exit(0)
