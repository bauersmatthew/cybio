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
from operator import xor

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
                        version='%(prog)s 1.1.0')
arg_parser.add_argument('isoforms',
                        type=valid_file,
                        help='The output file from cb-np_count-isoforms.py')
arg_parser.add_argument('filter',
                        type=str,
                        help=('The filter to be used. Format: \'A;B;C;...\' '
                              'where each component is \'[+-]feature\'. For '
                              'example: \'+11.4;-13.1\'. Basic (greedy) '
                              'boolean logic can also be used with operators '
                              '& (AND), | (OR), and ^ (XOR). For example: '
                              '\'+12.1^+12.12\'. The delimiter (;) is like an '
                              '& operator with parentheses surrounding its '
                              'arguments.'))
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

def tokenize(comp):
    """Tokenize one filter component."""
    l = []
    growing = ''
    for ch in comp:
        if ch in ('&', '|', '^'):
            l.append(growing)
            l.append(ch)
            growing = ''
        else:
            growing += ch
    l.append(growing)
    return tuple(l)

def validate(tokens):
    """Validate a tokenization.
    Raise exception on failure; exit silently on success."""
    for t in tokens[::2]:
        if len(t) < 2 or t[0] not in ('+', '-'):
            raise RuntimeError('Invalid filter token: \'{}\'!'.format(t))

args.filter = tuple(tokenize(c) for c in args.filter.split(args.delim))
for f in args.filter:
    validate(f)

def test_filter_component(iso, comp):
    """Test whether the given isoform matches the given filter component."""
    val = True
    op = '&'
    for i, tok in enumerate(comp):
        if i%2 == 1:
            op = tok
        else:
            sign = tok[0]
            region = tok[1:]
            if op == '&':
                if sign == '+': val &= (region in iso)
                else:           val &= (region not in iso)
            if op == '|':
                if sign == '+': val |= (region in iso)
                else:           val |= (region not in iso)
            if op == '^':
                if sign == '+': val ^= (region in iso)
                else:           val ^= (region not in iso)
    return val

def filter_one(line):
    line = line.rstrip()
    if not line:
        return ''
    isos = line.split('\t')[0].split(',')
    for f in args.filter:
        t = test_filter_component(isos, f)
        if not t:
            if args.negate:
                return line + '\n'
            else:
                return ''
    if args.negate:
        return ''
    else:
        return line + '\n'

# read isoform table
iso_counts = []
with open(args.isoforms) as fin:
    for line in fin:
        sys.stdout.write(filter_one(line))

sys.exit(0)
