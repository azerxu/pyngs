#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: ann.py
#
# annotation file formate
# **********************************************************************

class Ann(object):
    def __init__(self, mark, chrom, pos, ref, alt, refcount, altcount,
                 strand, gstart, gend, gene, exon):
        self.mark = mark
        self.chrom = chrom
        self.pos = int(pos) - 1
        self.ref = ref
        self.alt = alt
        self.refcount = int(refcount)
        self.altcount = int(altcount)
        self.bitstrand = 1 if strand == '+' else 2 # 1 as '+', 2 as '-'
        self.gstart = int(gstart) - 1
        self.gend = int(gend)
        self.gene = gene
        self.exon = exon

    @property
    def strand(self):
        return '+' if self.bitstrand == 1 else '-'

    def __repr__(self):
        return '\t'.join(
            map(str, (self.mark, self.chrom, self.pos+1, self.ref, self.alt,
                      self.refcount, self.altcount, self.strand, self.gstart+1,
                      self.gend, self.gene, self.exon)))


def parse(fname):
    with open(fname, 'r') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            if not line:
                continue

            items = line.split('\t')

            yield Ann(*items)


def read(fname):
    return parse(fname)

