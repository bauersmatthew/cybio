#!/usr/bin/env python

# Script: "cb-np_count-isoforms.py"
# Requires: Python 3
# Purpose: This script will count different transcript splicing isoforms present
#          in the given alignment.
# Author: Matthew Bauer

import argparse
import sys
import os.path

class BEDRecord:
    """One record from a BED file."""
    def __init__(self, line):
        fields = line.split('\t')
        if len(fields) < 3:
            raise RuntimeError('Invalid BED line!')
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        l = len(fields)
        self.nfields = l
        if l > 3: self.name = fields[3]
        if l > 4: self.score = fields[4]
        if l > 5: self.strand = fields[5]
        if l > 6: self.thickstart = fields[6]
        if l > 7: self.thickend = fields[7]
        if l > 8: self.color = fields[8]
        if l > 9:
            n_blocks = int(fields[9])
            block_sizes = [int(x) for x in fields[10][:n_blocks]]
            block_starts = [int(x) for x in fields[11][:n_blocks]]
            self.blocks = []
            for i in range(n_blocks):
                self.blocks.append((block_starts[i],
                                    block_starts[i]+block_sizes[i]))
        else:
            self.blocks = [(self.start, self.end)]

def tuple_overlap(t1, t2):
    if t2[1] <= t1[0] or t2[0] >= t1[1]:
        return 'none'
    elif t2[0] >= t1[0] and t2[1] <= t1[1]:
        return 'full'
    else:
        return 'partial'

def load_bed_file(path):
    """Load the BED file with the given path."""
    records = []
    header = ''
    with open(path) as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith('browser') or line.startswith('track'):
                header += '{}\n'.format(line)
            elif line:
                records.append(BEDRecord(line))
    return header, records

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

arg_parser = argparse.ArgumentParser(
    description=(
        'Count different transcript splicing isoforms present in an '
        'alignment.'),
    epilog='Results are written to STDOUT.')
arg_parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s 1.0.0')
arg_parser.add_argument('alignments',
                        help='A BED12 file of aligned nanopore reads.',
                        type=valid_file)
group = arg_parser.add_mutually_exclusive_group(required=True)
group.add_argument('-s', '--separated-annotation',
                   type=valid_file,
                   metavar='bed6_annotation',
                   help=('A BED6 file of gene structure annotations. '
                         'In this format, each exon is given its own BED '
                         'record, and exon names are taken from the '
                         'record name for each exon. Any exon blocks '
                         'present will be ignored.'))
group.add_argument('-b', '--blocked-annotation',
                   type=valid_file,
                   metavar='bed12_annotation',
                   help=('A BED12 file of gene structure annotations. '
                         'In this format, a gene should have only one '
                         'BED record associated with it; exons should be '
                         'described using the exon blocks fields of the '
                         'record.'))
args = arg_parser.parse_args()

_, aligned_reads = load_bed_file(args.alignments)
reference_path = args.separated_annotation
separated = True
if reference_path is None:
    reference_path = args.blocked_annotation
    separated = False
_, reference = load_bed_file(reference_path)

def init_count_structure():
    """The counting structure is different depending on the reference type."""
    global reference
    ret = {}
    for ref in reference:
        if separated:
            ret[ref.name] = 0
        else:
            ret[ref.name] = [0]*len(ref.blocks)
    return ret

# initialize count structure
counts = init_count_structure()

# go through each record
for read in aligned_reads:
    # count for this read
    my_count = init_count_structure()
    for block in read.blocks:
        for ref in reference:
            if separated and tuple_overlap(block, ref) != 'none':
                my_count[ref.name] = 1
            else:
                for rbi, refblock in enumerate(ref.blocks):
                    if tuple_overlap(block, refblock) != 'none':
                        my_count[ref.name][rbi] = 1
    # combine count for this read with total count
    for name in my_count:
        if separated:
            counts[name] += my_count[name]
        else:
            for i in range(len(my_count[name])):
                counts[name][i] += my_count[name][i]

# write results
for name in counts:
    if separated:
        sys.stdout.write('{}\t{}\n'.format(name, counts[name]))
    else:
        for exon in range(len(counts[name])):
            sys.stdout.write('{}:{}\t{}\n'.format(name, exon,
                                                  counts[name][exon]))
sys.exit(0)
