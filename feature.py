#!/usr/bin/env python


class Feature(object):
    """
    Minimalist class for handling feature entries from a GFF3 file.

    In GFF3, feature locations are described with 1-based indexing and closed
    intervals. The first 10 nucleotides of a sequence are represented using
    the notation [1-10].
    """
    def __init__(self, gff3_line):
        self._data = gff3_line
        fields = gff3_line.rstrip().split('\t')
        assert len(fields) == 9
        self.seqid = fields[0]
        self.ftype = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])

        self.attrs = dict()
        for keyvaluepair in fields[8].split(';'):
            key, value = keyvaluepair.split('=')
            self.attrs[key] = value

    def __repr__(self):
        return self._data.rstrip()

    def __len__(self):
        return self.end - self.start + 1

    def attr(self, key):
        if key not in self.attrs:
            return None
        return self.attrs[key]

    @property
    def label(self):
        return self.attr('Name')
