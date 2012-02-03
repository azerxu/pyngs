#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: splitpair.py
#
# merge paired fasta file or fastq file into one
# **********************************************************************

import sys


def split(fname, fmt='fq'):
    if fmt in ('fq', 'fastq'):
        from pyngs.biofile.fastq import parse
    elif fmt in ('fa', 'fna', 'fasta', 'fsa'):
        from pyngs.biofile.fasta import parse

    for idx, read in enumerate(parse(fname)):
        if idx % 2:                     # read2
            print >>sys.stderr, read
        else:                           # read1
            print >>sys.stdout, read


def main(args):
    if not args:
        print 'Usage: splitpair.py [-fq, -fa] fname1 fname2 ....'
        exit()

    fmt = 'fq'
    if args[0].startswith('-'):
        fmt = args[0][1:]
        args = args[1:]

    for arg in args:
        split(arg, fmt=fmt)


if __name__ == '__main__':
    main(sys.argv[1:])

