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
            block_starts = [b[0] for b in self.blocks]
            block_sizes = [b[1]-b[0] for b in self.blocks]
            ret += '\t{}\t{}\t{}'.format(
                            n_blocks,
                            ','.join(map(str, block_sizes)),
                            ','.join(map(str, block_starts)))
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

    def get_overlap(self, other):
        """Get the number of bases overlap between this record and another."""
        if self.start <= other.start and self.end > other.start:
            return self.end-other.start
        elif other.start <= self.start and other.end > self.start:
            return other.end-self.start
        return 0

    def get_coverage(self):
        """Get the total number of bases covered by all blocks."""
        return sum([b[1]-b[0] for b in self.blocks])

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

def get_most_overlapping(rec, test_recs):
    """Get the record in test_recs that rec overlaps with most."""
    most = max(test_recs, key=lambda tr: rec.get_overlap(tr))
    if rec.get_overlap(most) == 0:
        return None
    return most

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
tot_reads_killed = 0
tot_bases_begin  = 0
tot_bases_end    = 0
for read in annotated:
    tot_bases_begin += read.get_coverage()
    gene = get_most_overlapping(read, reference)
    if gene is None:
        tot_reads_killed += 1
        continue
    found_blocks = []
    for gb in gene.blocks:
        for rb in read.blocks:
            if ((rb[0] >= gb[0] and rb[0] < gb[1]) or
                (rb[1] > gb[0] and rb[1] <= gb[1])):
                found_blocks.append(gb)
                break
    if found_blocks:
        read.blocks     = found_blocks
        read.start      = gene.start+found_blocks[0][0]
        read.end        = gene.start+found_blocks[-1][1]
        read.thickstart = read.start
        read.thickend   = read.end
        sys.stdout.write('{}\n'.format(read))
        tot_bases_end += read.get_coverage()
    else:
        tot_reads_killed += 1

# maybe write statistics
if not args.quiet:
    n_total         = len(aligned_reads)
    n_annotated     = len(annotated)
    n_unannotated   = n_total-n_annotated
    n_bases_changed = tot_bases_end - tot_bases_begin
    say = sys.stderr.write
    say('{} reads processed.\n'.format(
        n_total))
    say('{} annotated ({}%).\n'.format(
        n_annotated,
        int(100*n_annotated/n_total)))
    say('{} unannotated ({}%).\n'.format(
        n_unannotated,
        int(100*n_unannotated/n_total)))
    qualifier = 'added' if n_bases_changed >= 0 else 'removed'
    say('{} bases {} overall ({}%; {} bases per read).\n'.format(
        abs(n_bases_changed),
        qualifier,
        int(100*abs(n_bases_changed)/tot_bases_begin),
        int(abs(n_bases_changed)/n_annotated)))
    say('{} annotated reads had no exons ({}%).\n'.format(
        tot_reads_killed,
        int(100*tot_reads_killed/n_annotated)))

sys.exit(0)
