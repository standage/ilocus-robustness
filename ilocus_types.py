#!/usr/bin/env python

from __future__ import print_function
import argparse
import re
import sys


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ids', metavar='FILE', type=argparse.FileType('r'),
                        help='if provided, only report iLocus_type for the '
                        'specified iLoci')
    parser.add_argument('locusgff3', type=argparse.FileType('r'),
                        help='iLocus GFF3 file')
    return parser


def main(args):
    loci2keep = None
    if args.ids:
        loci2keep = dict()
        for line in args.ids:
            locusid = line.strip()
            loci2keep[locusid] = True

    locus_types = dict()
    for record in args.locusgff3:
        if '\tlocus\t' not in record or 'iLocus_type' not in record:
            continue

        if loci2keep is not None:
            idmatch = re.search('Name=([^;\n]+)', record)
            assert idmatch, record
            locusid = idmatch.group(1)
            if locusid not in loci2keep:
                continue

        typematch = re.search('iLocus_type=([^;\n]+)', record)
        assert typematch, record
        locustype = typematch.group(1)
        if locustype not in locus_types:
            locus_types[locustype] = 0
        locus_types[locustype] += 1

    total = 0
    for locustype in locus_types:
        count = locus_types[locustype]
        print(locustype, count, sep=': ')
        total += count
    print('Total: ', total)


if __name__ == '__main__':
    main(get_parser().parse_args())
