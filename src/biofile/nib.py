#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: nib.py
# parse .nib format file
# **********************************************************************

# .nib files

# A .nib file describes a DNA sequence packing two bases into each byte.
# Each nib file contains only a single sequence. A nib file begins with a
# 32 bit signature which is 0x6BE93D3A in the architecture of the machine
# that created the file, and possibly a byte-swapped version of the same
# number on another machine. This is followed by a 32 bit number in the
# same format which describes the number of bases in the file. This is
# followed by the bases themselves packed two bases to the byte. The first
# base is packed in the high order 4 bits (nibble), the second base in the
# low order four bits. In C code:

#     byte = (base1<<4) + base2
#     The numerical values for the bases are:
#     0 - T,  1 - C,  2 - A,  3 - G,  4 - N (unknown)
#     The most significant bit in a nibble is set if the base is masked.

import os
import struct

# Constant Define
NIB_MAGIC = 0x6BE93D3A
NIB_MASK = 0b1111


# BASECODE = {0: 'T', 1: 'C', 2: 'A', 3: 'G', 4: 'N'}
BASECODE = 'TCAGN'


class Nib(object):
    def __init__(self, nibfile, handle, nbases, offset):
        self._filename = nibfile
        self._handle = handle
        self._nbases = nbases
        self._offset = offset

    @property
    def nbase(self):
        return self._nbases

    def __getitem__(self, idx):
        if isinstance(idx, int):
            if idx > self._nbases:
                return ''
            elif idx < 0:
                return ''
            else:
                loc = idx // 2
                # self._handle.seek(self._offset)
                num = self._handle.seek(self._offset + loc, 0).read(1)
                if idx % 2:
                    return get_base2(num)
                else:
                    return get_base1(num)
        elif isinstance(idx, tuple) or isinstance(idx, list):
            pass



def get_base1(num):
    base = num >> 4
    return BASECODE[base]


def get_base2(num):
    base = num & NIB_MASK
    return BASECODE[base]


def parse(nibfile):
    bases = []
    with open(nibfile, 'rb') as handle:
        signatue, bcount = struct.unpack('<ii', handle.read(8))
        assert signatue == NIB_MAGIC, '{0} is not a nib file'.format(nibfile)

        for i in xrange(bcount // 2):
            num, = struct.unpack('<B', handle.read(1))
            bases.append(get_base1(num))
            bases.append(get_base2(num))

        if not (bcount + 1) % 2:
            num, = struct.unpack('<B', handle.read(1))
            bases.append(get_base1(num))

    return ''.join(bases)


def read(nibfile):
    return parse(nibfile)

