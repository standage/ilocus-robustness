#!/usr/bin/env python


import re


class Match(object):
    """
    Very minimal class for handling records from vmatch output.

    This class is not an attempt at a robust representation of sequence
    alignments or robust handling of vmatch output. It is intentionally minimal
    to support a very limited subset of operations based on very specific
    assumptions about how vmatch was run.

    Important to note: vmatch reports alignment locations using 0-based
    indexing. The 0-based start position of each alignment is provided, along
    with the length of the alignment. The end position is not provided.
    """
    def __init__(self, line):
        line = line.strip()
        values = re.compile(' +').split(line)
        self.subject_length = int(values[0])
        self.subject_seqid = values[1]
        self.subject_pos = int(values[2])
        self.strand = values[3]
        self.query_length = int(values[4])
        self.query_seqid = values[5]
        self.query_pos = int(values[6])
        self.distance = int(values[7])
        self.evalue = float(values[8])
        self.score = int(values[9])
        self.identity = float(values[10])
        self.query = None

    @property
    def substart(self):
        """Start coordinate of the match on the subject."""
        return self.subject_pos

    @property
    def subend(self):
        """End coordinate of the match on the subject."""
        return self.subject_pos + self.subject_length

    def resolve(self, intervals):
        if intervals is None:
            return None

        valid_match = list()
        for interval in intervals:
            subject = interval.data
            overlap = (min(self.query.end, subject.end) -
                       max(self.query.start, subject.start))
            perc_query = float(overlap) / float(len(self.query))
            perc_subject = float(overlap) / float(len(subject))
            if perc_query >= 0.9 and perc_subject >= 0.9:
                valid_match.append(subject)

        if len(valid_match) == 0:
            return None
        return valid_match
