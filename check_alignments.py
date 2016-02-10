#!/usr/bin/env python


from __future__ import print_function
import argparse
import sys

import feature
import intervalforest
import match


def parse_iloci(infile):
    for line in infile:
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue

        feat = feature.Feature(line)
        if feat.ftype != 'locus':
            continue

        yield feat


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('q_gff3', type=argparse.FileType('r'),
                        help='iLocus annotation for query')
    parser.add_argument('s_gff3', type=argparse.FileType('r'),
                        help='iLocus annotation for subject')
    parser.add_argument('vmatch', type=argparse.FileType('r'),
                        help='vmatch output file')
    return parser


def main(args):
    # Store query iLoci by label
    print('Load query iLoci...', file=sys.stderr)
    q_iloci = dict()
    for locus in parse_iloci(args.q_gff3):
        q_iloci[locus.label] = locus

    # Store subject iLoci by label, and by interval for
    # efficient range-based queries
    print('Load subject iLoci...', file=sys.stderr)
    s_iloci = dict()
    forest = intervalforest.IntervalForest()
    for locus in parse_iloci(args.s_gff3):
        s_iloci[locus.label] = locus
        forest.add_feature(locus)

    print('Process vmatch matches...', file=sys.stderr)
    s_iloci_matched = dict()
    for line in args.vmatch:
        if line.startswith('#'):
            continue
        m = match.Match(line)
        m.query = q_iloci[m.query_seqid]
        overlapping_iloci = forest.search(m.subject_seqid, m.substart,
                                          m.subend)
        subject_iloci = m.resolve(overlapping_iloci)
        subject_label = None
        if subject_iloci is not None:
            for si in subject_iloci:
                s_iloci_matched[si.label] = m.query
            subject_label = ','.join([x.label for x in subject_iloci])
        print(m.query.label, subject_label, sep='\t')

    for interval in forest:
        locus = interval.data
        if locus.label not in s_iloci_matched:
            print(None, locus.label, sep='\t')


if __name__ == '__main__':
    main(args=get_parser().parse_args())
