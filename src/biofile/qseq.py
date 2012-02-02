#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: qseq.py
# **********************************************************************

"""Illumina QSEQ files processing
According to Illumina manual qseq files have the following format:

(1) Machine name: (hopefully) unique identifier of the sequencer.
(2) Run number: (hopefully) unique number to identify the run on the sequencer.
(3) Lane number: positive integer (currently 1-8).
(4) Tile number: positive integer.
(5) X: x coordinate of the spot. Integer (can be negative).
(6) Y: y coordinate of the spot. Integer (can be negative).
(7) Index: positive integer. No indexing should have a value of 1.
(8) Read Number: 1 for single reads; 1 or 2 for paired ends.
(9) Sequence
(10) Quality: the calibrated quality string.
(11) Filter: Did the read pass filtering? 0 - No, 1 - Yes.
"""

from fastq import Fastq
from xopen import xopen


PHRED64_FORMAT = set(("SOLEXA", "ILLUMINA", "PHRED64",
                      "ILLUMINA1.3", "ILLUMINA1.5",
                      "ILL1.3", "ILL1.5",
                      'X', 'I', 'J'))


def phred64to33():
    _64 = ''.join((chr(c) for c in xrange(64, 126)))
    _33 = ''.join((chr(c) for c in xrange(33, 33 + 62)))
    return maketrans(_64, _33)


def parse(qseqfile, fmt='I'):
    fmt = fmt.upper()
    handle = xopen(qseqfile, 'rb')
    table_64_to_33 = phred64to33()
    for line in handle:
        line = line.strip()
        if not line:
            continue

        (mach, runid, lane, tile, x, y, index ,readid,
         seq, qual, fil) = line.split('\t')

        # if fil value is 1 pass filter, 0 not
        fil = 'N' if fil == '1' else 'Y'

        if fmt in PHRED64:
            # trans phred64 quality to phred33 quality
            qual = qual.translate(table_64_to_33)

        name = '{0}:{1}:{2}:{3}:{4}:{5} {6}:{7}:{8}'.format(
            mach, runid, lane, tile, x, y, readid, fil, index)
        yield Fastq(name, seq, qual)


def read(qseqfile, fmt='S'):
    return parse(qseqfile, fmt=fmt)

