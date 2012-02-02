#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: trimer.py
#
# trim fastq file by quality score
# **********************************************************************
"""Trimer a tool to trimer NGS fastq format data by quality score

Usage: trimer.py [opts] fastqfile1 fastqfile2 fastqfile3 ....
       detail option as below:
       -m or --method str  trim algorithm: best or bwa [bwa]
       -t or --tag    str  output file name prefix, if not given then
                           use input file name as prefix default is ''
       -q or --qthres int  quality threshold value default is [20]
       -l or --lthres int  length threshold value default is [35]
       -f or --qtype  str  quality value type: S or I [S]
       -p or --pair        input files are paired
       -h or --help        show this help message
"""

import os
from pyngs.biofile.fastq import parse
from itertools import izip
import getopt


QTHRESHOLD = 20
LTHRESHOLD = 35
GOOD = 0                                # both reads pass trim
BAD1 = 1                                # read1 not pass trim
BAD2 = 2                                # read2 not pass trim
BOTH = 3                                # both reads are not pass trim


def _best(fastq, qthres=QTHRESHOLD):
    """All base quality must great than threshold

    Arguments:
    - `fastq`: Fastq Object
    - `qthres`: quality cutoff threshold
    """
    # maxlen record the max length value
    # best_start record start position for good
    # cur_start record the current start position
    maxlen = best_start = cur_start = 0

    for pos, qval in enumerate(fastq.qval):
        if qval < qthres:                # when qval not satisfy threshold
            if pos - cur_start > maxlen: # compare current length with _max
                best_start = cur_start   # record new max
                maxlen = pos - cur_start
            cur_start = pos + 1         # renew current start pos

    if len(fastq) - cur_start > maxlen: # check the last
        best_start = cur_start
        maxlen = len(fastq) - cur_start
    return best_start, maxlen           # fastq[best_start:best_start+maxlen]


def _bwa(fastq, qthres=QTHRESHOLD):
    """Calculate the max of continue base quality value great than threshold
    means: make sure max(sum(each - threshold))
    Arguments:
    - `fastq`: Fastq Object
    - `threshold`: cutoff threshold
    """
    # maxval record the max (sum(basequality - qthres)) value
    # best_start record start position for good
    # cur_start record the current start position
    # maxlen record the max base length after trim
    maxval = best_start = cur_start = maxlen = 0

    cumval = 0                          # record cum sum value
    for pos, qval in enumerate(fastq.qval):
        cumval += (qval - qthres)
        if cumval >= maxval:
            maxval = cumval             # record new max value
            best_start = cur_start
            maxlen = pos - cur_start + 1
        elif cumval < 0:
            cumval = 0
            cur_start = pos + 1
    return best_start, maxlen


def get_method(method='bwa'):
    """return a method to trim Fastq Object

    Arguments:
    - `method`:
    """
    assert method in ('bwa', 'best'), 'Unkown parameter:{0}'.format(method)
    return eval('_' + method)


def _trim_single(fname, method='bwa', qtype='S', qthres=QTHRESHOLD,
                lthres=LTHRESHOLD):
    """trim fastq by quality values

    Arguments:
    - `fname`: fastq file name
    - `method`: the way to trim Fastq object
    - `qtype`: the fastq qulity type
    - `qthres`: quality cutoff threshold
    - `lthres`: length cutoff threshold
    """
    _trim = get_method(method=method)
    for fq in parse(fname, fmt=qtype):
        start, length = _trim(fq, qthres)
        if length < lthres:
            yield BAD1, fq
        else:
            yield GOOD, fq[start:start+length]


def _trim_pair(pair1, pair2, method='bwa', qtype='S', qthres=QTHRESHOLD,
              lthres=LTHRESHOLD):
    _trim = get_method(method)
    fq1s = parse(pair1, fmt=qtype)
    fq2s = parse(pair2, fmt=qtype)
    for fq1, fq2 in izip(fq1s, fq2s):
        mark = 0
        start1, length1 = _trim(fq1, qthres)
        start2, length2 = _trim(fq2, qthres)
        if length1 < lthres:
            mark |= 1
        elif length2 < lthres:
            mark |= 2

        if mark == GOOD:
            yield GOOD, fq1[start1:start1+length1], fq2[start2:start2+length2]
        else:
            yield mark, fq1, fq2


def trim_single(fname, **kwargs):
    qthres = kwargs.get('qthres', QTHRESHOLD)
    lthres = kwargs.get('lthres', LTHRESHOLD)

    qtype = kwargs.get('qtype', 'S')
    method = kwargs.get('method', 'bwa')
    is_pair = kwargs.get('ispair', False)

    tag = kwargs.get('tag', '')

    if not tag:
        tag = os.path.splitext(os.path.basename(fname))[0]

    name_temp = '{0}.{1}.q{2}l{3}.{4}.fq'

    good = open(name_temp.format(tag, method, qthres, lthres, 'good'), 'w')
    bad = open(name_temp.format(tag, method, qthres, lthres, 'bad'), 'w')
    log = open('{0}.{1}.q{2}l{3}.log'.format(tag, method, qthres, lthres), 'w')

    for mark, fq in _trim_single(fname, method=method, qthres=qthres,
                                 qtype=qtype, lthres=lthres):
        if mark == GOOD:
            print >>good, fq
        else:
            print >>bad, fq
        print >>log, '\t'.join((str(mark), fq.name))

    for out in (good, bad, log):
        out.close()


def trim_pair(pair1, pair2, **kwargs):
    qthres = kwargs.get('qthres', QTHRESHOLD)
    lthres = kwargs.get('lthres', LTHRESHOLD)

    qtype = kwargs.get('qtype', 'S')
    method = kwargs.get('method', 'bwa')
    is_pair = kwargs.get('ispair', False)

    name_temp = '{0}.{1}.q{2}l{3}.{4}{5}.fq'

    tag = kwargs.get('tag', '')
    if not tag:
        tag = os.path.splitext(os.path.basename(pair1))[0]

    good1 = open(name_temp.format(tag, method, qthres, lthres, 'good', '1'), 'w')
    good2 = open(name_temp.format(tag, method, qthres, lthres, 'good', '2'), 'w')
    bad1 = open(name_temp.format(tag, method, qthres, lthres, 'bad', '1'), 'w')
    bad2 = open(name_temp.format(tag, method, qthres, lthres, 'bad', '2'), 'w')
    log = open('{0}.{1}.q{2}l{3}.log'.format(tag, method, qthres, lthres), 'w')

    for mark, fq1, fq2 in _trim_pair(pair1, pair2, method=method, qthres=qthres,
                                     qtype=qtype, lthres=lthres):
        if mark == GOOD:
            print >>good1, fq1
            print >>good2, fq2
        else:
            print >>bad1, fq1
            print >>bad2, fq2
        print >>log, '\t'.join((str(mark), fq1.name, fq2.name))

    for out in (good, bad, log):
        out.close()


def show_usage():
    print __doc__
    exit()


def main(argv):
    ispair = False

    kwargs = {
        'qtype': 'S',
        'tag': '',
        'qthres': QTHRESHOLD,
        'lthres': LTHRESHOLD,
        'method': 'bwa',
        'ispair': ispair
        }

    try:
        optlst, args = getopt.getopt(
            argv, 'hf:t:q:l:m:p', ['help', 'qtype', 'tag', 'qthres', 'lthres',
                                   'method', 'pair'])
        # parse arguments to trim data
        for opt, val in optlst:
            if opt in ('-h', '--help'): # show usage
                show_usage()
            elif opt in ('-t', '--tag'): # output prefix name
                kwargs['tag'] = val
            elif opt in('-m', '--method'): # trim method
                kwargs['method'] = val
            elif opt in ('-f', '--qtype'): # fastq quality format
                kwargs['qtype'] = val.upper()
            elif opt in ('-q', '--qthres'): # quality threshold
                kwargs['qthres'] = int(val)
            elif opt in ('-l', '--lthres'): # length threshold
                kwargs['lthres'] = int(val)
            elif opt in ('-p', '--pair'): # input file is paired or not
                ispair = True
                kwargs['ispair'] = ispair
    except getopt.GetoptError, e:
        show_usage()

    # no sff file given
    if not args:
        show_usage()

    if ispair:                          # process paired data
        for idx in xrange(0, len(args), 2):
            trim_pair(args[idx], args[idx+1], **kwargs)
    else:                               # process single data
        for arg in args:
            trim_single(arg, **kwargs)


# **********************************************************************
# main code block
# **********************************************************************
if __name__ == "__main__":
    import sys
    main(sys.argv[1:])

