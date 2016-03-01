#!/usr/bin/env python


from __future__ import print_function
import argparse
import intervaltree
import itertools
import sys

import feature
import match
import mapping


def resolve(match, intervals, logfile=None):
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
        overlap = (min(match.end, subject.end) -
                   max(match.start, subject.start))
        perc_query = float(overlap) / float(len(match.query))
        perc_subject = float(overlap) / float(len(subject))
        if logfile:
            qloc = '%s:%d-%d' % match.loc
            sloc = '%s:%d-%d' % subject.loc
            message = 'Candidate mapping: '
            message += 'query=%s[%s] ' % (match.query.label, qloc)
            message += 'subject=%s[%s] ' % (subject.label, sloc)
            message += 'qovlp=%.4lf sovlp=%.4lf' % (perc_query, perc_subject)
            print(message, file=logfile)
        if perc_query >= 0.9 and perc_subject >= 0.9:
            valid_mappings.append(subject)

    if len(valid_mappings) == 0:
        return None
    return valid_mappings


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
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout, metavar='OF',
                        help='output file; default is terminal (stdout)')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w'),
                        default=None, metavar='LF', help='print verbose '
                        'diagnostic messages to the specified file')
    parser.add_argument('q_gff3', type=argparse.FileType('r'),
                        help='iLocus annotation for query')
    parser.add_argument('s_gff3', type=argparse.FileType('r'),
                        help='iLocus annotation for subject')
    parser.add_argument('vmatch', type=argparse.FileType('r'), nargs='+',
                        help='vmatch output file(s)')
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
    q_iloci_matched = dict()
    locusmapping = mapping.Mapping()
    for line in itertools.chain(*args.vmatch):
        if line.startswith('#'):
            continue
        m = match.Match(line)
        m.query = q_iloci[m.query_seqid]
        seqid = m.subject_seqid
        if seqid not in intervals:
            continue
        overlapping_iloci = intervals[seqid].search(m.start, m.end)
        mapped_iloci = resolve(m, overlapping_iloci, args.logfile)
        locusmapping.add(m.query.label, mapped_iloci)
        q_iloci_matched[m.query.label] = True

    for locusid, mappedloci in locusmapping:
        print(locusid, mappedloci, sep='\t', file=args.outfile)

    for locusid in q_iloci:
        if locusid not in q_iloci_matched:
            print(locusid, None, sep='\t', file=args.outfile)

    for locusid in s_iloci:
        if locusid not in s_iloci_matched:
            print(None, locusid, sep='\t', file=args.outfile)


if __name__ == '__main__':
    main(args=get_parser().parse_args())
