#!/usr/bin/env python

from __future__ import print_function

class Mapping(object):
    def __init__(self):
        self.mapping = dict()

    def add(self, locusid, matches):
        if locusid not in self.mapping:
            self.mapping[locusid] = list()
        if matches is None:
            return
        for match in matches:
            self.mapping[locusid].append(match)

    def __iter__(self):
        for locusid in sorted(self.mapping):
            if len(self.mapping[locusid]) == 0:
                yield locusid, None
                continue
            yield locusid, [x.label for x in self.mapping[locusid]]

