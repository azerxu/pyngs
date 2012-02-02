#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: qseq2fq.py
#
# **********************************************************************

from pyngs.biofile.qseq import parse


def to(qseqfile, fmt='I'):
    for fq in parse(qseqfile, fmt=fmt):
        print fq


def main(args):
    fmt = 'S'
    if args and args[0].startswith('-'):
        fmt = args[0][1:]
        args = args[1:]

    if not args:
        print "Usage: qseq2fq.py [-format] qseqfile1 qseqfile2 ..."
        print "    format is qseq quality format, detail as be low:"
        print "    -S: Sanger (default)"
        print "    -I: Illumina"
        exit()

    map(lambda arg: to(arg, fmt), args)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

