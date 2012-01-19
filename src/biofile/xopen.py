#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: xopen.py
#
# Replacement for the "open" function that can also open
# files that have been compressed with gzip. If the filename ends with .gz,
# the file is opened with gzip.open(). If it doesn't, the regular open()
# is used. If the filename is '-', standard output (mode 'w') or input
# (mode 'r') is returned.
# **********************************************************************

import gzip
import sys

def xopen(fname, mode='r'):
    assert isinstance(fname, basestring)

    if fname == '-':
        if 'r' in mode:
            return sys.stdin
        elif 'w' in mode:
            return sys.stdout
        else:
            raise ValueError('Type wrong: {0}'.format(mode))

    if fname.endswith('.gz'):
        return gzip.open(fname, mode)

    return open(fname, mode)

