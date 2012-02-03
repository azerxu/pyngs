#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: agp.py
#
# process NCBI agp format file
# **********************************************************************

class Agp(object):
    def __init__(self, items):
        self.items = items

    @property
    def scaf(self):
        return self.items[0]

    @property
    def is_gap(self):
        return True if self.items[4] == 'N' else False

    @property
    def gaplen(self):
        return int(self.items[5])

    @property
    def scaf_start(self):
        return int(self.items[1]) - 1

    @property
    def scaf_end(self):
        return int(self.items[2])

    @property
    def pnum(self):
        return int(self.items[3])

    @property
    def type(self):
        return self.items[4]

    @property
    def contig(self):
        return self.items[5]

    @property
    def contig_start(self):
        return int(self.items[6]) - 1

    @property
    def contig_end(self):
        return int(self.items[7])

    @property
    def strand(self):
        return self.items[8]

    def __repr__(self):
        return '\t'.join(self.items + ([''] if self.is_gap else []))


def parse(agpfile):
    with open(agpfile) as handle:
        for line in handle:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            items = line.split('\t')
            agp = Agp(items)
            yield agp


def read(agpfile):
    return parse(agpfile)

