#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: mergepair.py
#
# merge paired fasta file or fastq file into one
# **********************************************************************

from itertools import izip

def merge(pair1, pair2, fmt='fq'):
    if fmt in ('fq', 'fastq'):
        from pyngs.biofile.fastq import parse
    elif fmt in ('fa', 'fna', 'fasta', 'fsa'):
        from pyngs.biofile.fasta import parse

    read1s = parse(pair1)
    read2s = parse(pair2)

    for read1, read2 in izip(read1s, read2s):
        print read1
        print read2


def main(args):
    if not args:
        print 'Usage: mergepair.py [-fq, -fa] pair1 pair2 ....'
        exit()

    fmt = 'fq'
    if args[0].startswith('-'):
        fmt = args[0][1:]
        args = args[1:]

    for i in xrange(0, len(args), 2):
        merge(args[i], args[i+1], fmt=fmt)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

