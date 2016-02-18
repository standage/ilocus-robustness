#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def parse_fasta(data):
    """Stolen shamelessly from http://stackoverflow.com/a/7655072/459780."""
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def format_seq(seq, linewidth=70, outstream=sys.stdout):
    """Print a sequence in a readable format."""
    if linewidth == 0 or len(seq) <= linewidth:
        print(seq, file=outstream)
        return

    i = 0
    while i < len(seq):
        print(seq[i:i+linewidth], file=outstream)
        i += linewidth


def get_parser():
    desc = 'Retrieve sequences by ID from Fasta data.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-o', '--out', type=argparse.FileType('w'),
                        default=sys.stdout, metavar='OF',
                        help='Output file; default is terminal (stdout)')
    parser.add_argument('-l', '--lwidth', type=int, default=70,
                        metavar='LW',
                        help='Max line width for sequences; default is 70 bp')
    parser.add_argument('idlist', type=argparse.FileType('r'))
    parser.add_argument('seqs', type=argparse.FileType('r'))
    return parser


def main(args):
    ids = dict()
    for line in args.idlist:
        seqid = line.rstrip()
        ids[seqid] = True

    for defline, seq in parse_fasta(args.seqs):
        seqid = defline[1:].split()[0]
        if seqid in ids:
            print(defline, file=args.out)
            format_seq(seq, linewidth=args.line_width, outstream=args.out)


if __name__ == '__main__':
    main(get_parser().parse_args())
