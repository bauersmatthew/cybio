#!/usr/bin/env python

# Script: "cb-np_draw-isoforms.py"
# Requires: Python 3
# Purpose: Generate a latex tikz drawing of all the isoforms found, including
#          their counts.
# Author: Matthew Bauer

import argparse
import sys
import os.path
import subprocess
import os

def valid_file(p):
    """Raise an exception if p is not a valid file; otherwise, return p."""
    if not os.path.isfile(p):
        raise RuntimeError("'{}' is not a valid input file!".format(p))
    return p

arg_parser = argparse.ArgumentParser(
    description=(
        'Generate a latex tikz drawing of all the isoforms found, including '
        'their counts.'),
    epilog='Results are written to STDOUT.')
arg_parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0.0')
arg_parser.add_argument('annotation',
                        help='A BED6 annotation of the gene of interest.',
                        type=valid_file)
arg_parser.add_argument('isoforms',
                        type=valid_file,
                        help='The output file from cb-np_count-isoforms.py.')
arg_parser.add_argument('-x', '--width',
                        type=float, required=False, default=13.0,
                        help=('The width of the image to be generated (cm). '
                              'Default: 13.0'))
arg_parser.add_argument('-y', '--height',
                        type=float, required=False, default=0.5,
                        help=('The height of each isoform (cm). '
                              'Default: 0.5'))
arg_parser.add_argument('-s', '--spacing',
                        type=float, required=False, default=0.25,
                        help=('The vertical space between each isoform (cm). '
                              'Default: 0.25'))
arg_parser.add_argument('-E', '--exon-style',
                        type=str, required=False,
                        default='blue!35!black,fill=blue!35!black')
arg_parser.add_argument('-I', '--intron-style',
                        type=str, required=False,
                        default='blue!35!black,ultra thick')
arg_parser.add_argument('-B', '--bold-style',
                        type=str, required=False,
                        default='red,fill=red')
arg_parser.add_argument('-b', '--bold',
                        type=str, required=False, default='',
                        help=('Comma-separated list of features to be '
                              'bolded.'))
arg_parser.add_argument('-t', '--template',
                        type=str, required=False, default='',
                        help=('Comma-separated list of features to be shown '
                              'in the template. Default: do not show a '
                              'template.'))
arg_parser.add_argument('-T', '--template-names',
                        type=str, required=False, default='',
                        help=('Names to be used for display in the '
                              'template. Default: use the internal feature '
                              'names. If given, must be a list the same '
                              'length as given with --template.'))
arg_parser.add_argument('-r', '--reverse',
                        action='store_true', default=False,
                        help=('Reverse the image (useful for showing the '
                              'antisense strand).'))
arg_parser.add_argument('-S', '--svg',
                        type=str, metavar='out.svg',
                        help=('Compile the output as a standalone TikZ '
                              'picture, producing an SVG file. Default: '
                              'only output the raw latex commands, starting '
                              'from \\begin{tikzpicture} and ending at '
                              '\\end{tikzpicture}.'))
args = arg_parser.parse_args()

args.bold = args.bold.split(',')
args.template = args.template.split(',')
tn = args.template_names.split(',')
if args.template != [''] and tn != ['']:
    args.template_names = {}
    for i in range(len(tn)):
        args.template_names[args.template[i]] = tn[i]

# read the annotation
features = {}
with open(args.annotation) as fin:
    for line in fin:
        line = line.strip()
        if not line:
            continue
        f = line.split('\t')
        features[f[3]] = (int(f[1]), int(f[2]))

# read in isoforms; while reading, find the range covered by all isos
min_start = 2**30 # very large; requires 31 bits
max_end = 0
isos = []
feats_present = set()
with open(args.isoforms) as fin:

    lines = list(fin.readlines())
    if args.template != ['']:
        lines.append('{}\t0'.format(','.join(args.template)))

    for line in lines:
        line = line.strip()
        if not line:
            continue
        f = line.split('\t')
        has = tuple(f[0].split(','))
        if(int(f[1]) != 0):
            isos.append((has, int(f[1])))
        for n in has:
            feats_present.add(n)
            min_start = min(min_start, features[n][0])
            max_end = max(max_end, features[n][1])

# calculate tikz coordinates
midpoint = min_start + ((max_end-min_start)//2)
scaling = args.width/(max_end-min_start)
draw_coords = {}
for f in feats_present:
    c = features[f]
    if not args.reverse:
        draw_coords[f] = (scaling*(c[0]-midpoint), scaling*(c[1]-midpoint))
    else:
        draw_coords[f] = (
            -scaling*(c[1]-midpoint),
            -scaling*(c[0]-midpoint))

# begin drawing
w = lambda s: sys.stdout.write(s+'\n')
fout = None
if args.svg is not None:
    fout = open(args.svg+'.tex', 'w')
def w(s):
    if fout is not None:
        fout.write(s+'\n')
    else:
        sys.stdout.write(s+'\n')
if fout is not None:
    w('\\documentclass[tikz,convert={{outfile={}}}]{{standalone}}'.format(
        args.svg))
    w('\\begin{document}')
w('\\begin{tikzpicture}')
s = lambda l: sorted(l, key=lambda x: draw_coords[x][0])
def draw_exon(name, y):
    style = args.exon_style
    if name in args.bold: style = args.bold_style
    w('\\draw[{}] ({},{}) rectangle ({},{});'.format(
        style,
        draw_coords[name][0], y+(args.height/2),
        draw_coords[name][1], y-(args.height/2)))
def draw_intron(e1n, e2n, y):
    w('\\draw[{}] ({},{}) -- ({},{});'.format(
        args.intron_style,
        draw_coords[e1n][1], y,
        draw_coords[e2n][0], y))
def draw(iso, y):
    f = s(iso[0])
    for i in range(len(f)-1):
        draw_exon(f[i], y)
        draw_intron(f[i], f[i+1], y)
    last = f[-1]
    draw_exon(last, y)
    w('\\node[anchor=west, right] at ({},{}) {{{}}};'.format(
        draw_coords[last][1] + (args.width/10), y,
        iso[1]))

# draw template
y = 0.0
if args.template != ['']:
    args.template = s(args.template)
    for i in range(len(args.template)-1):
        name = args.template[i]
        name2 = args.template[i+1]
        draw_exon(name, y)
        draw_intron(name, name2, y)
        c_s = draw_coords[name][0]
        c_e = draw_coords[name][1]
        draw_name = name
        if args.template_names:
            draw_name = args.template_names[name]
        if i%2 == 0:
            w('\\node[anchor=south, above] at ({},{}) {{{}}};'.format(
                c_s+((c_e-c_s)/2), y+args.height,
                draw_name))
        else:
            w('\\node[anchor=north, below] at ({},{}) {{{}}};'.format(
                c_s+((c_e-c_s)/2), y-args.height,
                draw_name))

    name = args.template[-1]
    draw_exon(name, y)
    c_s = draw_coords[name][0]
    c_e = draw_coords[name][1]
    draw_name = name
    if args.template_names:
        draw_name = args.template_names[name]
    if (i+1)%2 == 0:
        w('\\node[anchor=south, above] at ({},{}) {{{}}};'.format(
            c_s+((c_e-c_s)/2), y+args.height,
            draw_name))
    else:
        w('\\node[anchor=north, below] at ({},{}) {{{}}};'.format(
            c_s+((c_e-c_s)/2), y-args.height,
            draw_name))

    y -= 3*(args.height + args.spacing)

# draw each isoform
for iso in isos:
    draw(iso, y)
    y -= args.height + args.spacing
w('\\end{tikzpicture}')

if fout is not None:
    w('\\end{document}')
    fout.close()
    subprocess.run(['pdflatex', '--shell-escape', '{}.tex'.format(args.svg)])
    os.remove('{}.tex'.format(args.svg))
    os.remove('{}.aux'.format(args.svg))
    os.remove('{}.log'.format(args.svg))
    os.remove('{}.pdf'.format(args.svg))
