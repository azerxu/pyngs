#!/usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: fastq.py
# **********************************************************************

'''Parse Fastq File library

Reference: SURVEY AND SUMMARY The Sanger FASTQ file format for sequences
           with quality scores, and the Solexa/Illumina FASTQ variants
           Peter J. A. Cock, Christopher J. Fields, Naohisa Goto,
           Michael L. Heuer and Peter M. Rice
Published online 16 December 2009
FASTQ DEFINITION
@title and optional description
sequence line(s)
+optional repeat of title line
quality line(s)

example:
@SRR014849.1 EIXKN4201CFU84 length=93
GGGGGGGGGGGGGGGGCTTTTTTTGTTTGGAACCGAAAGG
GTTTTGAATTTCAAACCCTTTTCGGTTTCCAACCTTCCAA
AGCAATGCCAATA
+SRR014849.1 EIXKN4201CFU84 length=93
3+&$#"""""""""""7F@71,’";C?,B;?6B;:EA1EA
1EA5’9B:?:#9EA0D@2EA5’:>5?:%A;A8A;?9B;D@
/=<?7=9<2A8==
'''

# *********************************************************************
# Description, OBF name    ASCII characters    Quality score
#                          Range     Offset    Type     Range
# ------------------------------------------------------------
# Sanger standard
#   fastq-sanger          33–126     33        PHRED   0 to 93
# Solexa/early Illumina
#   fastq-solexa          59–126     64        Solexa -5 to 62
# Illumina 1.3+
#   fastq-illumina        64–126     64        PHRED   0 to 62
# ------------------------------------------------------------

# *********************************************************************
# Illumina sequence identifiers
# Sequences from the Illumina software use a systematic identifier:
# @HWUSI-EAS100R:6:73:941:1973#0/1
# HWUSI-EAS100R	the unique instrument name
# 6	flowcell lane
# 73	tile number within the flowcell lane
# 941	'x'-coordinate of the cluster within the tile
# 1973	'y'-coordinate of the cluster within the tile
# #0	index number for a multiplexed sample (0 for no indexing)
# /1	the member of a pair, /1 or /2 (paired-end or mate-pair reads only)
# *********************************************************************

# fastq-sanger, fastq-illumina Qual define
# Qphred = -10 * log10(Pe)
# Pe = 10 ** (Qphred/-10)

# fastq-solexa Qual define
# Qsolexa = -10 * log10(Pe / (1-Pe))
# Pe = 1 / (10 ** (Qsolexa/10) + 1)

# *****************************************************************************
#  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
#  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
#  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
#  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
#  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
#  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
#  |                         |    |        |                              |                     |
# 33                        59   64       73                            104                   126
# *****************************************************************************

# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
#     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
#     (Note: See discussion above).

# With Casava 1.8 the format of the '@' line has changed:
# @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
#
# EAS139	the unique instrument name
# 136	the run id
# FC706VJ	the flowcell id
# 2	flowcell lane
# 2104	tile number within the flowcell lane
# 15343	'x'-coordinate of the cluster within the tile
# 197393	'y'-coordinate of the cluster within the tile
# 1	the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
# Y	Y if the read fails filter (read is bad), N otherwise
# 18	0 when none of the control bits are on, otherwise it is an even number
# ATCACG	index sequence


from string import maketrans
from xopen import xopen
import gzip
import sys

# Quality OFFSET
PHRED33_OFFSET = 33                       # Standard offset
PHRED64_OFFSET = 64                       # Solexa or Illumina offset


# fastq format set
PHRED33_TYPE = set(('S', 'SANGER', 'PHRED33'))
PHRED64_TYPE = set(("SOLEXA", "ILLUMINA", "PHRED64",
                      "ILLUMINA1.3", "ILLUMINA1.5",
                      "ILL1.3", "ILL1.5",
                      'X', 'I', 'J'))


class Fastq(object):
    def __init__(self, name='', seq='', qual=''):
        """Fastq format
        name: read name and other annotation
        seq: sequence
        qual: quality values and the quality should be phred33 quality
        """
        self.name = name
        self.seq = seq
        self.qual = qual

    def __len__(self):
        """calc Fastq sequence length"""
        return len(self.seq)

    def __getitem__(self, idx):
        """Used for fastq slice"""
        return Fastq(name=self.name, seq=self.seq[idx],
                     qual=self.qual[idx])

    def __repr__(self):
        """for print and output"""
        return '@{0}\n{1}\n+\n{2}'.format(self.name, self.seq, self.qual)

    @property
    def qval(self):
        """Get quality number quality"""
        return tuple((ord(qval) - PHRED33_OFFSET) for qval in self.qual)

    @property
    def pval(self):
        """Get p value"""
        # fastq-sanger, fastq-illumina Qual define
        # Qphred = -10 * log10(Pe)
        # Pe = 10 ** (Qphred/-10)
        return (pow(10, -0.1 * q) for q in self.qval)


# deco function to generator trans
def dtrans(func):
    table = func()
    def _trans(qual):
        """translate qulaity"""
        return qual.translate(table)
    return _trans


@dtrans
def phred64to33():
    _64 = ''.join((chr(c) for c in xrange(64, 126)))
    _33 = ''.join((chr(c) for c in xrange(33, 33 + 62)))
    return maketrans(_64, _33)


@dtrans
def phred33to64():
    _64 = ''.join((chr(c) for c in xrange(64, 126)))
    _33 = ''.join((chr(c) for c in xrange(33, 33 + 62)))
    return maketrans(_33, _64)


def parse(fname, qtype='S'):
    """parse fastq file and return a iterator
    standard is a mark to show whether format to trans to standard
    """
    seq = ''
    qual = ''
    name = ''
    slen = qlen = 0

    if qtype in PHRED64_TYPE:
        need_trans = True
        trans = phred64to33
    else:
        need_trans = False

    is_seq_block = False                # True as Seq block, False Qual block

    handle = xopen(fname, 'r')
    # read head lines to check is or not fastq file
    for line in handle:
        line = line.rstrip()
        if not line or line.startswith('#'):
            continue
        if not line.startswith('@'):
            raise ValueError('{0} is not in fastq format'.format(fname))
        break                           # quit cycle

    # check is a empty fastq file or not
    if not line:
        return

    is_seq_block = True
    name = line[1:]

    for line in handle:
        line = line.rstrip()            # trim right endof \n \r
        if not line:                    # ignore blank line
            continue
        if is_seq_block:                # deal with seq block
            if line.startswith('+'):    # next is qual block
                is_seq_block = False
            else:                       # deal with seq block
                seq += line
                slen += len(line)
        else:                           # deal with quality block
            if qlen > slen:             # check qual length <= seq length
                raise ValueError('Error while Parsing {0}'.format(name))

            if line.startswith('@'):    # switch to sequence block
                # at beginning of next fastq
                if seq and slen == qlen:
                    if need_trans:
                        qual = trans(qual)
                    yield Fastq(name, seq, qual)
                    seq = ''
                    qual = ''
                    name = line[1:]
                    is_seq_block = True    # next is seq block
                elif not seq and not qual: # start to generate fastq
                    name = line[1:]
                    is_seq_block = True    # next is seq block
                else:                      # just a qual line begin with @
                    qual += line           # renew quality value
                    qlen += len(line)
            else:
                qual += line
                qlen += len(line)

    # yield last fastq record
    if name or seq:
        if slen != qlen:                    # check the last fastq record
            raise ValueError('parsing wrong with {0}'.format(name))

        if need_trans:                  # trans qual
            qual = trans(qual)

        yield Fastq(name, seq, qual)


def read(fname, qtype='S'):
    """read a fastq record from fastq file"""
    try:
        return parse(fname, qtype=qtype).next()
    except StopIteration:
        raise ValueError("Empty fastq file: {0}".format(fname))


if __name__ == '__main__':
    import sys
    for arg in sys.argv[1:]:
        for fq in parse(arg):
            print fq

