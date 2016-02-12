#!/usr/bin/env python


class Feature(object):
    """
    Very minimal class for handling feature entries from a GFF3 file.

    This class is not an attempt at a robust representation of genome features
    or robust handling of GFF3 data. It is intentially minimal to support a
    very limited subset of operations.

    Important to note: in GFF3, feature locations are described with 1-based
    indexing and closed intervals. For example, the first 10 nucleotides of a
    sequence are represented using the notation [1-10]. Internal to this
    implementation, however, locations are stored as 0-based half-open
    intervals, where the notation [0, 10) would be used to refer to the same
    10 nucleotides.
    """
    def __init__(self, gff3_line):
        self._data = gff3_line
        fields = gff3_line.rstrip().split('\t')
        assert len(fields) == 9
        self.seqid = fields[0]
        self.ftype = fields[2]
        self.start = int(fields[3]) - 1
        self.end = int(fields[4])

        self.attrs = dict()
        for keyvaluepair in fields[8].split(';'):
            key, value = keyvaluepair.split('=')
            self.attrs[key] = value

    def __repr__(self):
        return self._data.rstrip()

    def __len__(self):
        return self.end - self.start

    def attr(self, key):
        if key not in self.attrs:
            return None
        return self.attrs[key]

    @property
    def label(self):
        return self.attr('Name')

    @property
    def loc(self):
        return (self.seqid, self.start + 1, self.end)
