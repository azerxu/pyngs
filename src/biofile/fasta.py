#!/usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: fasta.py
# Parsing fasta format file
# **********************************************************************

"""Module contain two method (parse and read) and a Fasta Class:
function parse: return a Fasta object iterator
function read: return the first fasta record
class Fasta: easy way to deal with fasta record
"""

from xopen import xopen                 # get read gzip file support
LINE_WIDTH = 60                         # each line contain bases


class Fasta(object):
    def __init__(self, name='', seq=''):
        self.name = name
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        """seq in multi-line"""
        seqs = ['>{0}'.format(self.name)]
        for i in xrange(0, len(self.seq), LINE_WIDTH):
            seqs.append(self.seq[i:i+LINE_WIDTH])
        return '\n'.join(seqs)

    def __str__(self):
        """seq in one line"""
        return '>{0}\n{1}'.format(self.name, self.seq)

    def __getitem__(self, idx):
        """used for fasta slice"""
        return Fasta(name=self.name, seq=self.seq[idx])


def parse(fname):
    """Parse multi fasta records file and return a Fasta Object iterator"""
    name = ''
    seq = []
    handle = xopen(fname, 'r')
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name or seq:
                yield Fasta(name, ''.join(seq))
            name = line[1:]
            seq = []
        else:
            seq.append(line)

    if name or seq:
        yield Fasta(name, ''.join(seq))


def read(fname):
    """read fasta record from file"""
    try:
        return parse(fname).next()
    except StopIteration:
        raise ValueError, 'Fasta file: {0} is Empty'.format(fname)

