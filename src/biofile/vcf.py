#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: libvcf.py
# reference: http://www.1000genomes.org/wiki/Analysis/
# Variant%20Call%20Format/vcf-variant-call-format-version-41
# **********************************************************************
"""Read and Parse VCF (Variant Call Format) version 4.1"""

import re

# File meta-information is included after the ## string, often as
# key=value pairs.
# A single 'fileformat' field is always required, must be the first
# line in the file

# INFO fields should be described as follows (all keys are required):
# ##INFO=<ID=ID,Number=number,Type=type,Description=”description”>
RE_INFO = re.compile('##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<desc>[^"]+)">')

class Info(object):
    def __init__(self, id, number, type, desc):
        self._id = id
        self._number = number
        self._type = type
        self._desc = desc

    def __repr__(self):
        return '##INFO=<ID={0},Number={1},Type={2},Description="{3}">'\
               .format(self._id, self._number, self._type, self._desc)


# FILTERs that have been applied to the data should be described as follows:
# ##FILTER=<ID=ID,Description=”description”>
RE_FILTER = re.compile('##FILTER=<ID=[^,]+,Description="[^"]+">')

class Filter(object):
    def __init__(self, id, desc):
        self._id = id
        self._desc = desc

    def __repr__(self):
        return '##FILTER=<ID={0},Description="{1}">'.format(self._id, self._desc)

# Likewise, Genotype fields specified in the FORMAT field should be described
# as below:
# ##FORMAT=<ID=ID,Number=number,Type=type,Description=”description”>
RE_FORMAT = re.compile('##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<desc>[^"]+)">')

class Format(object):
    def __init__(self, id, number, type, desc):
        self._id = id
        self._number = number
        self._type = type
        self._desc = desc

    def __repr__(self):
        return '##FORMAT=<ID={0},Number={1},Type={2},Description="{3}">'\
               .format(self._id, self._number, self._type, self._desc)

# The header line names the 8 fixed, mandatory columns. These columns are as follows:
# #CHROM
# POS
# ID
# REF
# ALT
# QUAL
# FILTER
# INFO

class VcfRecord(object):
    def __init__(self, source, chrom, pos, id_, ref, alt, qual, filter_,
                 info, **kwargs):
        self.source = source
        self.chrom = chrom
        self.pos = pos
        self.id = id_
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_
        self.__dict__.update(**kwargs)
        self.deal_info(info)

    def deal_info(self, info):
        items = info.split(';')
        self.is_indel = False
        if items[0] == 'INDEL':
            self.is_indel = True
            items = items[1:]
        self.__dict__.update(**dict([item.split('=') for item in items if item]))

    @property
    def genotype(self):
        if 'PL' not in self.__dict__:
            return '{0}/{1}'.format(self.ref, self.ref)
        bases = [self.ref] + self.alt.split(',')
        types = []
        bn = len(bases)                 # bn: bases length
        for i in xrange(bn):
            for j in xrange(i, bn):
                types.append((bases[i], bases[j]))
        pls = self.PL.split(',')
        pls = [eval(each) for each in pls]
        _min = min(pls)
        idx  = pls.index(_min)
        return '{0}/{1}'.format(*types[idx])

    def __repr__(self):
        return '<VcfRecord>'


class Header(object):
    def __init__(self, items):
        self.items = items


class VCF(object):
    """Variant Call Format"""
    def __init__(self, handle, start, samples, **kwargs):
        self._handle = handle
        self._start = start
        self._samples = samples

    def reset(self):
        self._handle.seek(self._start)

    def __iter__(self):
        return self

    def next(self):
        for line in self._handle:
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line:
                continue
            try:
                return parse_vcfrecord(line, self._samples)
            except:
                raise StopIteration
        raise StopIteration


def _parse_meta(line, (regex, kls)):
    match = regex.match(line)
    if not match:
        return None
    dic = match.groupdict()
    return kls(**dic)


def parse_vcfrecord(line, samples):
    items = line.split('\t')
    chrom, pos, id_, ref, alt, qual, filter_, info = items[:8]
    _format = items[8].split(':')
    for i, item in enumerate(items[9:]):
        source = samples[i]
        data = item.split(':')
        format_ = dict(zip(_format, data))
        yield VcfRecord(source=source, chrom=chrom, pos=pos, id_=id_,
                        ref=ref, alt=alt, qual=qual, filter_=filter_,
                        info=info, **format_)


def parse(fname):
    handle = open(fname, 'r')
    for line in handle:                 # deal with head info block
        line = line.strip()
        if not line:
            continue
        if not line.startswith('##'):   #
            break

    if line.startswith('#'):        # header
        samples = line.split('\t')[9:]
    else:
        raise ValueError, 'No head lines'

    start = handle.tell()

    return VCF(handle, start, samples)


def read(fname):
    return parse(fname)

