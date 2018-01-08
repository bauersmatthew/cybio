#!/usr/bin/env python

# Script: "cb-np_count-isoforms.py"
# Requires: Python 3
# Purpose: This script, given the output of cb-np_count-isoforms.py, will
#          calculate the pairwise correlations for each pair of features.
# Author: Matthew Bauer

import argparse
import sys
import os.path
import re
from scipy import stats

def get_all_matches(l, expr):
    """Get all full REGEX matches to expr in l."""
    ret = []
    for e in l:
        m = re.match(expr, e)
        if m is not None and m.start() == 0 and m.end() == len(e):
            ret.append(e)
    return ret

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

arg_parser = argparse.ArgumentParser(
    description=('Calculate pairwise correlations for each pair of features.'),
    epilog='Results are written to STDOUT.')
arg_parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0.0')
arg_parser.add_argument('isoforms',
                        help='The output file from cb-np_count-isoforms.py.',
                        type=valid_file)
arg_parser.add_argument('-t', '--test',
                        action='append', type=str, metavar='test_pair',
                        help=('A description of pairs to calculate '
                              'correlations for. Format: first|second. '
                              'The first and second strings are REGEX '
                              'expressions. This argument can be passed more '
                              'than once to describe multiple test pairs.'))
args = arg_parser.parse_args()

# read isoform table; construct set of all features
iso_counts = []
features = set()
with open(args.isoforms) as fin:
    for line in fin:
        line = line.rstrip()
        if line:
            fields = line.split('\t')
            f = fields[0].split(',')
            c = int(fields[1])
            iso_counts.append(f, c)
            for feat in f:
                features.add(feat)

# collect features of interest
tests = []
for t in arg_parser.test:
    sides = t.split('|')
    if len(sides) != 2:
        sys.stderr.write('Improper test: \'{}\''.format(t))
        sys.exit(-1)
    s0_matches = get_all_matches(features, sides[0])
    s1_matches = get_all_matches(features, sides[1])
    for s0m in s0_matches:
        for s1m in s1_matches:
            tests.append((s0m, s1m))

def get_isos_with_and_without(f, counts):
    """Get isoform count lines with and without the given feature."""
    isos_with = []
    isos_without = []
    for c in counts:
        if f in c[0]:
            isos_with.append(c)
        else:
            isos_without.append(c)
    return isos_with, isos_without

def do_test(t):
    """Conduct one test."""
    yes0, no0 = get_isos_with_and_without(t[0], iso_counts)
    yes0yes1, yes0no1 = get_isos_with_and_without(t[1], yes0)
    no0yes1, no0no1 = get_isos_with_and_without(t[1], no0)

    # conduct chi-squared
    #    | y0 n0 |
    # ---+-------+
    # y1 | yy ny |
    # n1 | yn nn |
    # ---+-------+
    # (look at wiki)

    # observed values
    o_yy, o_yn, o_ny, o_nn = len(yes0yes1), len(yes0no1), len(no0yes1), \
        len(no0no1)
    # row/column totals
    tot = o_yy+o_yn+o_ny+o_nn
    r1t = o_yy+o_ny
    r2t = o_yn+o_nn
    c1t = o_yy+o_yn
    c2t = o_ny+o_nn
    # expected values
    e_yy = c1t*r1t/total
    e_ny = c2t*r1t/total
    e_yn = c1t*r2t/total
    e_nn = c2t*r2t/total
    # chi-squared partials
    x2_yy = pow(o_yy-e_yy, 2)/e_yy
    x2_ny = pow(o_ny-e_ny, 2)/e_ny
    x2_yn = pow(o_yn-e_yn, 2)/e_yn
    x2_nn = pow(o_nn-e_nn, 2)/e_nn
    x2 = x2_yy+x2_ny+x2_yn+x2_nn
    # compute p-value
    p = 1-stats.chi2.cdf(x2, 1)
    return p

# conduct tests, print results
sys.stdout.write('first\tsecond\tp-value\n')
for t in test:
    p = do_test(t)
    sys.stdout.write('{}\t{}\t{}\n'.format(t[0], t[1], do_test(t))
    sys.stdout.write('{}\t{}\t{}\n'.format(t[1], t[0], do_test((t[1], t[0])))
