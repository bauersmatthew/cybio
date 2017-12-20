#!/usr/bin/env python

# Script: "cb-np_smooth-to-annotation.py"
# Requires: Python 3
# Purpose: Due to the high error rate of nanopore sequencing, aligned reads
#          often map imperfectly to the gene annotation. To make processing 
#          these alignments slightly easier, this script will try to fix obvious
#          mistakes in the reads, such as tiny gaps or insertions.
# Author: Matthew Bauer

import argparse
import sys
import copy
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
            block_sizes = [int(x) for x in fields[10].split(',')[:n_blocks]]
            block_starts = [int(x) for x in fields[11].split(',')[:n_blocks]]
            self.blocks = []
            for i in range(n_blocks):
                self.blocks.append((block_starts[i],
                                    block_starts[i]+block_sizes[i]))
        else:
            self.blocks = [(self.start, self.end)]

    def str(self):
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
            block_starts = [b[0] for b in self.blocks]
            block_sizes = [b[1]-b[0] for b in self.blocks]
            ret += '\t{}\t{}\t{}'.format(n_blocks,
                                         ','.join(block_sizes),
                                         ','.join(block_starts))
        return ret

    def check_overlap(self, other):
        """Check how/if the given record overlaps this one.

        Return 'none' if they do not overlap at all;
               'partial' if other extends outside self;
               'full' if self fully eclipses other."""
        if other.end <= self.start or other.start >= self.end:
            return 'none'
        elif other.start >= self.start and other.end <= self.end:
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

class BinaryRegionalMask:
    """A mask that is either "on" or "off" over an integer region."""
    def __init__(self, on=[]):
        self.on = copy.copy(on)
        self._sort()

    def get_on_intervals(self):
        return copy.copy(self.on)

    def get_coverage(self):
        return sum(r[1]-r[0] for r in self.on)

    def get_inverse(self):
        ret = copy.copy(self)
        ret.invert()
        return ret

    def invert(self):
        """Invert the mask over the region containing all "on" intervals."""
        new_on = []
        for i in range(len(self.on)-1):
            r1 = self.on[i]
            r2 = self.on[i+1]
            new_on.append((r1[1], r2[0]))
        self.on = new_on

    def grow_seeds(self, other):
        """Grow all regions in this BRM that are subsets of regions in the other
        BRM to be the same size of those regions of which they are subsets."""
        to_add = []
        for reg in other.on:
            for seed in self.on:
                if seed[0] >= reg[0] and seed[1] <= reg[1]:
                    to_add.append(reg)
                    break
        for reg in to_add:
            self += reg

    def __add__(self, other):
        """Add two masks together (like bitwise or).
        Also works for mask and tuple."""
        cpy = copy.copy(self)
        if isinstance(other, tuple):
            if len(other) != 2:
                raise ValueError()
            cpy.on.append(other)
            cpy._sort()
            while cpy._combine_right(): pass
        else: # assume is a BinaryRegionalMask-like object
            for r in other.on:
                cpy += r
        return cpy

    def __sub__(self, other):
        """Subtract one mask from another (like bitwise xor).
        Also works for mask and tuple."""
        cpy = copy.copy(self)
        if isinstance(other, tuple):
            if len(other) != 2:
                raise ValueError()
            new_on = []
            for r in cpy.on:
                if other[0] >= r[0] and other[0] < r[1]:
                    # keep section on the left
                    new_on.append((r[0], other[0]))
                if other[1] > r[0] and other[1] <= r[1]:
                    # keep section on the right
                    new_on.append((other[1], r[1]))
            cpy.on = new_on
            cpy._remove_artifacts()
        else: # assume is a BinaryRegionalMask-like object
            for r in other.on:
                cpy -= r
        return cpy, n_removed

    def _combine_right(self):
        """Combine regions when necessary from left-to right ONE time.
    
        Assumes in correct order by starting position.
        Return True if anything was changed; false otherwise."""
        changed = False
        for i in range(len(self.on)-1):
            new_on = []
            r1 = self.on[i]
            r2 = self.on[i+1]
            if r1[1] > r2[0]: # NOT EQUAL TO!
                new_on.append((r1[0], max(r1[1], r2[1])))
                changed = True
            else:
                new_on.append(r1)
                new_on.append(r2)
        self.on = new_on
        return changed

    def _sort(self):
        """Sort in increasing order by starting position."""
        self.on.sort(key=lambda x: x[0])

    def _remove_artifacts(self):
        """Remove empty ranges of the form (a, a)."""
        self.on = [r for r in self.on if r[0] != r[1]]

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

arg_parser = argparse.ArgumentParser(
    description=(
        'Due to the high error rate of nanopore sequencing, aligned reads '
        'often map imperfectly to the gene annotation. To make processing '
        'these alignments slightly easier, this script will try to fix obvious '
        'mistakes in the reads, such as tiny gaps or insertions.'),
    epilog="The new 'smoothed' alignment are written to STDOUT.")
arg_parser.add_argument('alignments',
                        help='A BED12 file of aligned nanopore reads.',
                        type=valid_file)
arg_parser.add_argument('annotations',
                        help='A BED12 file of gene structure annotations.',
                        type=valid_file)
arg_parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s 1.0.0')
arg_parser.add_argument('-r', '--remove-unannotated',
                        action='store_true',
                        help=('Throw out reads that are not accounted for in '
                              'the given annotation. Default: do nothing.'))
arg_parser.add_argument('-S', '--process-spillovers',
                        action='store_true',
                        help=('Process reads that extend outside the annotated '
                              'gene regions normally. Default: treat them as '
                              'unannotated.'))
arg_parser.add_argument('-q', '--quiet',
                        action='store_true',
                        help=('Do not output statistics. Default: write '
                              'various interesting statistics to STDERR.'))
args = arg_parser.parse_args()

alignment_header, aligned_reads = load_bed_file(args.alignments)
_, reference = load_bed_file(args.annotations)

# check the reference for validity (nothing should overlap)
for i1 in range(len(reference)):
    for i2 in range(i1+1, len(reference)):
        if reference[i1].check_overlap(reference[i2]) != 'none':
            sys.stderr.write('Invalid annotation: records must not overlap!\n')
            sys.exit(-1)

# echo the header from the alignment file
sys.stdout.write(alignment_header)

# separate annotated from unannotated alignments
annotated = []
for read in aligned_reads:
    matched = False
    for gene in reference:
        overlap = gene.check_overlap(read)
        if ((overlap == 'partial' and args.process_spillovers) or
            overlap == 'full'):
            annotated.append(read)
            matched = True
            break
    if (not matched) and (not args.remove_unannotated):
        # echo immediately
        sys.stdout.write('{}\n'.format(read))

# process each annotated record; write output when done with each
tot_removed = 0
tot_removed_ratio = 0
tot_added = 0
tot_added_ratio = 0
tot_netchange = 0
for read in annotated:
    read_brm = BinaryRegionalMask(read.blocks)
    for gene in reference:
        gene_brm = BinaryRegionalMask(gene.blocks)
        # remove insertions in intron regions
        cov_before = read_brm.get_coverage()
        read_brm -= gene_brm.get_inverse()
        cov_after = read_brm.get_coverage()
        tot_removed += -(cov_after-cov_before)
        tot_removed_ratio += -(cov_after-cov_before)/cov_before
        tot_netchange += (cov_after-cov_before)
        # fill in gaps in exon regions
        cov_before = cov_after
        read_brm.grow_seeds(gene_brm)
        cov_after = read_brm.get_coverage()
        tot_added += (cov_after-cov_before)
        tot_added_ratio += (cov_after-cov_before)/cov_before
        tot_netchange += (cov_after-cov_before)
    read.blocks = read_brm.get_on_intervals()
    sys.stdout.write('{}\n'.format(read))

# maybe write statistics
if not args.quiet:
    n_total = len(aligned_reads)
    n_annotated = len(annotated)
    n_unannotated = n_total - n_annotated
    sys.stderr.write('Totals...\n')
    sys.stderr.write('\t{} annotated ({}%)\n'.format(
        n_annotated, 100*n_annotated//n_total))
    sys.stderr.write('\t{} unannotated ({}%)\n'.format(
        n_unannotated, 100*n_unannotated//n_total))
    sys.stderr.write('Processing averages...\n')
    sys.stderr.write('\t{:.1f} bases removed ({:.1}% of read)\n'.format(
        tot_removed/n_annotated, 100*tot_removed_ratio/n_annotated))
    sys.stderr.write('\t{:.1f} bases added ({:.1}% of read)\n'.format(
        tot_added/n_annotated, 100*tot_added_ratio/n_annotated))
    sys.stderr.write('\t({:+.1f} bases net)\n'.format(
        tot_netchange/n_annotated))

sys.exit(0)
