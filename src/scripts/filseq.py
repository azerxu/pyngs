#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: filseq.py
#
# filter sequnce by three step:
# 1. remove contain N reads
# 2. remove average quality below quality threshold when reads in fastq format
# 3. remove contain adapter containment
# 4. get pass filter reads
# **********************************************************************

from pyngs.lib import calign
from pyngs.util import revcom
from itertools import izip


FADAPTER1 = 'GATCGGAAGAGCGGTTTTCAGCAGGAATGCCGAG'
FADAPTER2 = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'

RADAPTER1 = revcom(FADAPTER1)
RADAPTER2 = revcom(FADAPTER2)

FASTQ_FMT = ('fq', 'fastq')
FASTA_FMT = ('fa', 'fasta', 'fsa', 'fna')


def find_best_alignment(adapters, seq, max_error_rate=0.1, overlap=15):
    best_result = None
    best_matches = 0
    best_index = None
    for index, adapter in enumerate(adapters):
        # try find exact match first
        pos = seq.find(adapter)
        if pos >= 0:                    # find exact result
            result = (0, len(adapter), pos, pos + len(adapter), len(adapter), 0)
            return result, index

        # use sw_algor to find result
        result = calign.globalalign_locate(adapter, seq, max_error_rate)

        (astart, astop, rstart, rstop, matches, errors) = result
        if astart == 0 and astop == len(adapter): # full adapter matched
            return result, index

        length = astop - astart

        if length < overlap:
            continue

        if matches > best_matches:
            best_result = result
            best_matches = matches
            best_index = index
    return (best_result, best_index)


def filpair(tag, pair1, pair2, ffmt='fastq', qfmt='S', overlap=10, qthred=20):
    """filpair: filter paired-end data
    tag: output tag name
    pair1: paired-end-1 filename
    pair2: paired-end-2 filename
    ffmt: file format type default is 'fastq'
    qfmt: read quality format, default is 'S' means Sanger format
    qthred: used for quality filter, default is 20
    overlap: used for filter adapter, default is 10
    """
    adps = (FADAPTER1, FADAPTER2, RADAPTER1, RADAPTER2)
    if ffmt == 'qseq':
        from pyngs.biofile.qseq import parse
    elif ffmt in FASTQ_FMT:
        from pyngs.biofile.fastq import parse
    elif ffmt in FASTA_FMT:
        from pyngs.biofile.fasta import parse
    else:
        raise ValueError('Unkown file format: {0}'.format(ffmt))

    is_fq = False

    if ffmt in FASTA_FMT:
        read1s = parse(pair1)
        read2s = parse(pair2)
        ext = '.fa'
    else:
        read1s = parse(pair1, fmt=qfmt)
        read2s = parse(pair2, fmt=qfmt)
        ext = '.fq'
        is_fq = True

    # pass filter
    pass1 = open(tag + '.pass1' + ext, 'w')
    pass2 = open(tag + '.pass2' + ext, 'w')

    # reads contain N
    n1 = open(tag + '.n1' + ext, 'w')
    n2 = open(tag + '.n2' + ext, 'w')

    # reads quality is lower than quality threshold
    if is_fq:
        low1 = open(tag + '.low1' + ext, 'w')
        low2 = open(tag + '.low2' + ext, 'w')

    # reads contain adpater containment
    adap1 = open(tag + '.adap1' + ext, 'w')
    adap2 = open(tag + '.adap2' + ext, 'w')

    # record adapter filter recorder detail
    adplog = open(tag + '.adp.log', 'w')

    # process step:
    # 1. remove contain N reads
    # 2. remove average quality below quality threshold when reads in fastq format
    # 3. remove contain adapter containment
    # 4. get pass filter reads
    for read1, read2 in izip(read1s, read2s):
        # process reads contain N
        if 'N' in read1.seq or 'N' in read2.seq:
            print >>n1, read1
            print >>n2, read2
            continue

        # process reads average quality score less thand quality threshold
        if is_fq and (sum(read1.qval) / len(read1) < qthred or sum(read2.qval) / len(read2) < qthred):
            print >>low1, read1
            print >>low2, read2
            continue

        # process reads contation adapter containment
        res_read1, idx_read1 = find_best_alignment(adps, read1.seq, overlap=overlap)
        res_read2, idx_read2 = find_best_alignment(adps, read2.seq, overlap=overlap)
        if res_read1 or res_read2:      # at least one read contain adapter
            print >>adap1, read1
            print >>adap2, read2

            if res_read1:                 # read1 contain adapter
                (astart, astop, rstart, rstop, matches, errors) = res_read1
                length = rstop - rstart
                adapter = adps[idx_read1]
                print >>adplog, '\t'.join(map(str, (
                    '1', read1.name, idx_read1, adapter, read1.seq[rstart:rstop],
                    length, astart+1, astop, rstart+1, rstop, matches)))

            if res_read2:               # read2 contain adapter
                (astart, astop, rstart, rstop, matches, errors) = res_read2
                length = rstop - rstart
                adapter = adps[idx_read2]
                print >>adplog, '\t'.join(map(str, (
                    '2', read2.name, idx_read2, adapter, read2.seq[rstart:rstop],
                    length, astart+1, astop, rstart+1, rstop, matches)))
            continue

        # if reads come to here then reads means pass the filter
        print >>pass1, read1
        print >>pass2, read2

    # close output handle
    for out in (pass1, pass2, adap1, adap2, adplog):
        out.close()

    if is_fq:                           # close low quality output handle
        for out in (low1, low2):
            out.close()
    # filter all done here!


def main(args):
    tag = 'fil'
    if args and args[0].startswith('-'):
        tag = args[0][1:]
        args = args[1:]

    if len(args) < 2:
        print 'Usage: filseq.py [-tag] pair1, pair2'
        exit()

    ffmt = 'qseq'
    qfmt = 'I'
    for i in xrange(0, len(args), 2):
        filpair(tag, args[i], args[i+1], ffmt=ffmt, qfmt=qfmt)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

