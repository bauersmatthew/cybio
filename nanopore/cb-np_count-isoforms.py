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
        if l > 4: self.score = int(fields[4])
        if l > 5: self.strand = fields[5]
        if l > 6: self.thickstart = int(fields[6])
        if l > 7: self.thickend = int(fields[7])
        if l > 8: self.color = fields[8]
        if l > 9:
            n_blocks = int(fields[9])
            block_sizes = list(map(int, fields[10].split(',')[:n_blocks]))
            block_starts = list(map(int, fields[11].split(',')[:n_blocks]))
            self.blocks = []
            for i in range(n_blocks):
                self.blocks.append((self.start+block_starts[i],
                                    self.start+block_starts[i]+block_sizes[i]))

    def __str__(self):
        """Store this record as a single BED line (no newline)."""
        ret = '{}\t{}\t{}'.format(self.chrom, self.start, self.end)
        l = self.nfields
        if l > 3: ret += '\t{}'.format(self.name)
        if l > 4: ret += '\t{}'.format(self.score)
        if l > 5: ret += '\t{}'.format(self.strand)
        if l > 6: ret += '\t{}'.format(self.thickstart)
        if l > 7: ret += '\t{}'.format(self.thickend)
        if l > 8: ret += '\t{}'.format(self.color)
        if l > 9:
            n_blocks = len(self.blocks)
            block_starts = [b[0]-self.start for b in self.blocks]
            block_sizes = [b[1]-b[0] for b in self.blocks]
            ret += '\t{}\t{}\t{}'.format(
                            n_blocks,
                            ','.join(map(str, block_sizes)),
                            ','.join(map(str, block_starts)))
        return ret

def get_overlap_interval(a, b):
    """Get the number of bases overlap between two intervals."""
    return max(0, min(a[1], b[1])-max(a[0], b[0]))

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

features = {} # name -> interval
for gene in reference:
    if not separated:
        for i, block in enumerate(gene.blocks):
            name = '{}:{}'.format(gene.name, i+1)
            if name in features:
                raise RuntimeError('Feature names must be unique!')
            features[name] = block
    else:
        if gene.name in features:
            raise RuntimeError('Feature names must be unique!')
        features[gene.name] = (gene.start, gene.end)

tallies = {} # sorted tuple of names -> count
for read in aligned_reads:
    found = []
    for f_name in features:
        for r_block in read.blocks:
            if get_overlap_interval(r_block, features[f_name]) > 0:
                if f_name not in found:
                    found.append(f_name)
    if found:
        tup = tuple(sorted(found))
        if tup not in tallies:
            tallies[tup] = 0
        tallies[tup] += 1

# write results
for pair in sorted(tallies.items(), key=lambda i: i[1], reverse=True):
    tup = pair[0]
    num = pair[1]
    sys.stdout.write('{}\t{}\n'.format(','.join(tup), num))

sys.exit(0)
