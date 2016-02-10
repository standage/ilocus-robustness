#!/usr/bin/env python

import intervaltree


class IntervalForest(object):
    """
    A set of interval trees, one for each sequence.

    This `interval forest` is a very simple data structure for storing
    intervals annotated on multiple genomic sequences. The interval forest
    stores a separate interval tree for each sequence. Accordingly, overlap and
    containment queries (using the `search` method) must specify not only a
    range but also a sequence identifier.

    Note: ranges and intervals in Python are implemented using half-open
    notation. That is, [0, 10) refers to the first 10 elements of a list,
    sequence, etc. This is consistent with the notation used by the BED format
    but is inconsistent with GTF/GFF3 which used a 1-based index and closed
    interval notation: in other words, [1, 10] to refer to the same interval.
    """

    def __init__(self):
        self.trees = dict()

    def add_feature(self, feature, fix_notation=True):
        """
        Add a feature to the interval forest.

        Feature is expected to have the following attributes:
        - seqid (string)
        - start (integer)
        - end (integer)

        It is assumed the feature is using GFF3-style interval notation, which
        is adjusted to be consistent with Python and the implementation of the
        interval tree. Set `fix_notation=False` if the provided feature is
        already using half-open interval notation.
        """
        if feature.seqid not in self.trees:
            self.trees[feature.seqid] = intervaltree.IntervalTree()
        start = feature.start
        if fix_notation:
            start -= 1
        end = feature.end
        self.trees[feature.seqid].addi(start, end, feature)

    def search(self, seqid, start, end, strict=False):
        """
        Search for features overlapping the specified interval.

        When `strict=True`, a containment query (rather than an overlap query)
        is performed.
        """
        if seqid not in self.trees:
            return None
        return self.trees[seqid].search(start - 1, end, strict=strict)

    def __iter__(self):
        """Iterate through all features in the interval forest."""
        for seqid in self.trees:
            for locus in self.trees[seqid]:
                yield locus
