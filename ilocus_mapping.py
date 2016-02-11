#!/usr/bin/env python


from __future__ import print_function
import argparse
import intervaltree
import sys

import feature
import match


def resolve(match, intervals):
    """
    Compute iLocus mapping.

    The `match` object represents the alignment of an iLocus from assembly A to
    assembly B, and `match.query` is a reference to the corresponding iLocus
    annotation record. The `intervals` object is a set of intervals from
    assembly B that overlap with the alignment. Each interval object contains a
    reference to an iLocus annotation record from assembly B.

    We designate a mapping of the iLocus from assembly A to an iLocus from
    assembly B *valid* if the alignment constitutes at least 90% reciprocal
    overlap of both iLoci.
    """
    if intervals is None:
        return None

    valid_mappings = list()
    for interval in intervals:
        # Here, `subject` is a candidate target of the iLocus query.
        subject = interval.data
        overlap = (min(match.query.end, subject.end) -
                   max(match.query.start, subject.start))
        perc_query = float(overlap) / float(len(self.query))
        perc_subject = float(overlap) / float(len(subject))
        if perc_query >= 0.9 and perc_subject >= 0.9:
            valid_mappings.append(subject)

    if len(valid_mappings) == 0:
        return None
    return valid_match


def parse_iloci(infile):
    """Very simple GFF3 parser."""
    for line in infile:
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue

        feat = feature.Feature(line)
        if feat.ftype != 'locus':
            continue

        yield feat


def get_parser():
    """Specify a command-line interface for this script."""
    parser = argparse.ArgumentParser()
    parser.add_argument('q_gff3', type=argparse.FileType('r'),
                        help='iLocus annotation for query')
    parser.add_argument('s_gff3', type=argparse.FileType('r'),
                        help='iLocus annotation for subject')
    parser.add_argument('vmatch', type=argparse.FileType('r'),
                        help='vmatch output file')
    return parser


def main(args):
    """Define the main method of this script."""

    # Store query iLoci by label
    print('Load query iLoci...', file=sys.stderr)
    q_iloci = dict()
    for locus in parse_iloci(args.q_gff3):
        q_iloci[locus.label] = locus

    # Store subject iLoci by label, and by interval for
    # efficient range-based queries
    print('Load subject iLoci...', file=sys.stderr)
    s_iloci = dict()
    intervals = dict()
    for locus in parse_iloci(args.s_gff3):
        s_iloci[locus.label] = locus
        if locus.seqid not in intervals:
            intervals[locus.seqid] = intervaltree.IntervalTree()
        intervals[locus.seqid].addi(locus.start, locus.end, locus)

    print('Process vmatch alignments...', file=sys.stderr)
    s_iloci_matched = dict()
    for line in args.vmatch:
        if line.startswith('#'):
            continue
        m = match.Match(line)
        m.query = q_iloci[m.query_seqid]
        seqid = m.subject_seqid
        if seqid not in intervals:
            print(m.query.label, None, sep='\t')
            continue
        overlapping_iloci = intervals[seqid].search(m.substart, m.subend)
        mapped_iloci = m.resolve(overlapping_iloci)
        subject_label = None
        if mapped_iloci is not None:
            for ilocus in mapped_iloci:
                s_iloci_matched[ilocus.label] = m.query
            subject_label = ','.join([x.label for x in mapped_iloci])
        print(m.query.label, subject_label, sep='\t')

    for seqid in intervals:
        for interval in intervals[seqid]:
            locus = interval.data
            if locus.label not in s_iloci_matched:
                print(None, locus.label, sep='\t')


if __name__ == '__main__':
    main(args=get_parser().parse_args())