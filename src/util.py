#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: util.py
# common useful tools collections
# **********************************************************************

import sys
import time
import os
from string import maketrans


TRANS_TABLE = maketrans('ACGTRYMKWSNacgtrymkwsn', 'TGCAYRKMWSNtgcayrkmwsn')


def get_basename(fname):
    return os.path.splitext(os.path.basename(fname))[0]


def reverse(seq):
    return ''.join(reversed(seq))


def complement(seq):
    return ''.join(seq.translate(TRANS_TABLE))


def revcom(seq):
    return ''.join(reversed(seq.translate(TRANS_TABLE)))


# **********************************************************************
# show process with time stamp
# **********************************************************************
def deco_show_info(info):
    sinfo = 'Starting {0} ...'.format(info)
    einfo = '{0} finished!'.format(info)
    def _show(func):
        def __run(*args, **kwargs):
            print >>sys.stderr, time.strftime('%Y-%m-%d %H:%M:%S'), sinfo
            res = func(*args, **kwargs)
            print >>sys.stderr, time.strftime('%Y-%m-%d %H:%M:%S'), einfo
            return res
        return __run
    return _show


def show_info(info):
    """Show info with time date display
    Arguments:
    - `info`: the text to output
    """
    print >>sys.stderr, time.strftime('%Y-%m-%d %H:%M:%S'), info

