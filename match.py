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
    def start(self):
        """Start coordinate of the match on the subject."""
        return self.subject_pos

    @property
    def end(self):
        """End coordinate of the match on the subject."""
        return self.subject_pos + self.subject_length

    @property
    def loc(self):
        return (self.subject_seqid, self.start + 1, self.end)
