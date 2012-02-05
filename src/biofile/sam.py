#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: filsam.py
# **********************************************************************
"""The SAM Format Specification
SAM stands for Sequence Alignment/Map format. It is a TAB-delimited text
format consisting of a header section, which is optional, and an alignment
section. If present, the header must be prior to the alignments. Header
lines start with ‘@’, while alignment lines do not. Each alignment line has
11 mandatory fields for essential alignment information such as mapping
position, and variable number of optional fields for flexible or aligner
specific information.
"""

# import re


# Constant variable
TAB = '\t'

FORWARD_STRAND = '+'     # forward strand
REVERSE_STRAND = '-'     # reverse strand

MASK_MULTI_FRAG = 0x1    # template having multiple fragments in sequencing
MASK_PROP_ALIGN = 0x2    # each fragment properly aligned according to the
                         # aligner
MASK_UNMAPPED = 0x4      # fragment unmapped
MASK_MATE_UNMAPPED = 0x8 # next fragment in the template unmapped
MASK_REVERSE = 0x10      # SEQ being reverse complemented
MASK_MATE_REVERSE = 0x20 # SEQ of the next fragment in the template being
                         # reversed
MASK_READ1 = 0x40        # the first fragment in the template
MASK_READ2 = 0x80        # the last fragment in the template
MASK_SECONDARY = 0x100   # secondary alignment
MASK_FAILQC = 0x200      # not passing quality controls
MASK_DUPLICATE = 0x400   # PCR or optical duplicate


# CIGAR_RE = '\*|(?:[0-9]+[MIDNSHPX=])+'


# class Cigar(object):
#     """
#     CIGAR string: \*|([0-9]+[MIDNSHPX=])+
#     The CIGAR operations are given in the following table (set ‘*’ if
#     unavailable):
#     ------------------------------------------------------------------
#     Op   BAM   Description
#     ------------------------------------------------------------------
#     M     0    alignment match (can be a sequence match or mismatch)
#     I     1    insertion to the reference
#     D     2    deletion from the reference
#     N     3    skipped region from the reference
#     S     4    soft clipping (clipped sequences present in SEQ)
#     H     5    hard clipping (clipped sequences NOT present in SEQ)
#     P     6    padding (silent deletion from padded refence)
#     =     7    sequence match
#     X     8    sequence mismatch
#     ------------------------------------------------------------------
#     * H can only be present as the first and/or last operation
#     * S may only have H operations between them and the ends of the cigar
#       string
#     * For mRNA-to-genome alignment, an N operation represents an intron.
#       For othe types of alignments, the interpretation of N is not defined
#     * Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
#     """
#     def __init__(self, cigar):
#         assert re.match(CIGAR_RE, cigar), 'cigar string must be form by {0}'\
#                .format(CIGAR_RE)
#         self._cigar = cigar

#     def __repr__(self):
#         return self._cigar

#     def __len__(self):
#         flen = 0                        # full length
#         clen = 0                        # current length
#         for c in self._cigar:
#             if c == '*':
#                 return 0
#             elif c in 'MIS=X':
#                 flen += clen
#                 clen = 0
#             elif c.isdigit():
#                 clen = clen * 10 + int(c)
#             elif c in 'DNSHP':              # not count
#                 clen = 0
#         return flen


class Flag(object):
    """
    ======================================================================
    FLAG: bitwise FLAG. Each bit is explained in the following table:
    ====    ==============================================================
    Bit     Description
    ====    ==============================================================
    0x1     template having multiple fragments in sequencing
    0x2     each fragment properly aligned according to the aligner
    0x4     fragment unmapped
    0x8     next fragment in the template unmapped
    0x10    SEQ being reverse complemented
    0x20    SEQ of the next fragment in the template being reversed
    0x40    the first fragment in the template
    0x80    the last fragment in the template
    0x100   secondary alignment
    0x200   not passing quality controls
    0x400   PCR or optical duplicate
    ======================================================================
    """
    def __init__(self, flag):
        self._flag = int(flag)

    def __repr__(self):
        return str(self._flag)

    @property
    def qstrand(self):                  # query reads strand
        """0x10: SEQ being reverse complemented"""
        return REVERSE_STRAND if self.flag & MASK_REVERSE else FORWARD_STRAND

    @property
    def mstrand(self):                  # mate reads strand
        """0x20: SEQ of the next fragment in the template being reversed"""
        return (REVERSE_STRAND if self.flag & MASK_MATE_REVERSE
                else FORWARD_STRAND)


class Sam(object):
    """
    ===========================================================================
    Col Field Type    Regexp/Range             Brief description
    ===========================================================================
     1  QNAME String  [!-?A-~]{1,255}          Query template NAME
     2  FLAG  Int     [0,216 -1]               bitwise FLAG
     3  RNAME String  \*|[!-()+-<>-~][!-~]*    Reference sequence NAME
     4  POS   Int     [0,229 -1]               1-based leftmost mapping
                                               POSition
     5  MAPQ  Int     [0,28 -1]                MAPping Quality
     6  CIGAR String  \*|([0-9]+[MIDNSHPX=])+  CIGAR string
     7  RNEXT String  \*|=|[!-()+-<>-~][!-~]*  Ref. name of the mate/next
                                               fragment
     8  PNEXT Int     [0,2^29 -1]              Position of the mate/next
                                               fragment
     9  TLEN  Int     [-2^29 +1,2^29 -1]       observed Template LENgth
    10  SEQ   String  \*|[A-Za-z=.]+           fragment SEQuence
    11  QUAL  String  [!-~]+                   ASCII of Phred-scaled base
                                               QUALity+33
    ===========================================================================
    """
    def __init__(self, qname, flag, rname, pos, mapq, cigar,
                 rnext, pnext, tlen, seq, qual, *tags):
        self.qname = qname              # Query template Name
        self.flag = int(flag)           # bitwise Flag
        self.rname = rname              # Reference Sequence Name
        # 1-based leftmost mapping Position trans to 0-based
        self.pos = int(pos) - 1 if pos != '*' else -1
        self.mapq = int(mapq)           # MAPping Quality
        self.cigar = cigar              # cigar string
        self.rnext = rnext              # Ref. name of the mate/next fragment
        # Position of the mate/next fragment
        self.pnext = int(pnext) -1 if pnext != '*' else -1
        self.tlen = int(tlen)           # observed Template LENgth
        self.seq = seq                  # fragment Sequence
        self.qual = qual                # ascii of phred-scaled base quality+33
        self._tags = {}                 # store tags of sam mapping
        self._parse_tags(tags)

    def __repr__(self):
        pos = '*' if self.pos == -1 else self.pos + 1
        pnex = '*' if self.pnext == -1 else self.pnext + 1
        return '\t'.join(
            map(str, [self.qname, self.flag, self.rname, pos,
                      self.mapq, self.cigar, self.rnext, pnex,
                      self.tlen, self.seq, self.qual] + self.tags))

    def _parse_tags(self, tags):
        for tag in tags:
            key, tp, val = tag.split(':')
            self._tags[key] = (tp, val)

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]

        key = key.upper()
        if key in self._tags:
            tp, val = self._tags[key]
            if tp == 'i':
                val = int(val)
            return val
        return None

    def get_cigar(self):
        l = 0

        for c in self.cigar:
            if c.isdigit():
                l = l*10 + int(c)
            else:
                yield c, l
                l = 0

    @property
    def qstrand(self):                  # query reads strand
        """0x10: SEQ being reverse complemented"""
        return REVERSE_STRAND if self.flag & MASK_REVERSE else FORWARD_STRAND

    @property
    def mstrand(self):                  # mate reads strand
        """0x20: SEQ of the next fragment in the template being reversed"""
        return (REVERSE_STRAND if self.flag & MASK_MATE_REVERSE
                else FORWARD_STRAND)

    @property
    def aend(self):
        """aligned end position of the read on the reference genome.
        Return None if not available"""
        pass

    @property
    def alen(self):
        """aligned length of read on the reference genome.
        Return None if not available"""
        pass

    @property
    def is_paired(self):
        """0x1: template having multiple fragments in sequencing"""
        return self.flag & MASK_MULTI_FRAG

    @property
    def is_read1(self):
        """0x40: the first fragment in the template"""
        return self.flag & MASK_READ1

    @property
    def is_read2(self):
        """0x80: the last fragment in the template"""
        return self.flag & MASK_READ2

    @property
    def is_reverse(self):
        """0x10: SEQ being reverse complemented"""
        return self.flag & MASK_REVERSE

    @property
    def is_secondary(self):
        """True if not primary alignment
        0x100: secondary alignment
        Bit 0x100 marks the alignment not to be used in certain analyses
        when the tools in use are aware of this bit
        """
        return self.flag & MASK_SECONDARY

    @property
    def is_proper_pair(self):
        """0x2: each fragment properly aligned according to the aligner
        True if read is paired in sequence"""
        return self.flag & MASK_PROP_ALIGN

    @property
    def is_unmapped(self):
        """0x4: fragment unmapped
        Bit 0x4 is the only reliable place to tell whether the fragment is
        unmapped. If 0x4 is set, no assumptions can be made about RNAME,
        POS, CIGAR, MAPQ, bits 0x2, 0x10 and 0x100 and the bit 0x20 of the
        next fragment in the template.
        """
        return self.flag & MASK_UNMAPPED

    @property
    def mate_is_reverse(self):
        """True is read is mapped to reverse strand"""
        return self.flag & MASK_MATE_REVERSE

    @property
    def mate_is_unmapped(self):
        """True if the mate is unmapped
        0x8: next fragment in the template unmapped"""
        return self.flag & MASK_MATE_UNMAPPED

    @property
    def is_failqc(self):
        """0x200: not passing quality controls"""
        return self.flag & MASK_FAILQC

    @property
    def is_duplicate(self):
        """0x400: PCR or optical duplicate"""
        return self.flag & MASK_DUPLICATE

    def qstart(self):
        """start index of the aligned query portion of the sequence
        (0-based, inclusive)"""
        pass

    def qqual(self):
        """aligned query sequence quality values (None if not present)"""
        pass

    def qend(self):
        """end index of the aligned query portion of the sequence
        (0-based, exclusive)"""
        pass

    def qlen(self):
        """length of the aligned query sequence"""
        pass

    def query(self):
        """aligned portion of the read and excludes any flanking bases that
        were soft clipped (None if not present"""
        pass

    @property
    def rlen(self):
        """length of the read. Returns 0 if not given"""
        return len(self.seq) if self.seq else 0

    @property
    def readid(self):
        """get read id: 1 if read is first fragment, 2 the last one"""
        if self.flag & MASK_READ1:
            return 1
        if self.flag & MASK_READ2:
            return 2
        return 1

    @property
    def tags(self):
        _tags = []
        for key in self._tags:
            tp, val = self._tags.get(key)
            _tags.append('{0}:{1}:{2}'.format(key, tp, val))
        return _tags


class SamFile(object):
    def __init__(self, filename, header, handle, offset):
        self.filename = filename
        self.header = header
        self._offset = offset
        self._handle = handle

    def __iter__(self):
        return self

    def reset(self):
        self._handle(self._offset, 0)

    def next(self):
        while True:
            line = self._handle.readline()
            if not line:
                raise StopIteration
            line = line.rstrip()
            if line:
                return _parse_line(line)

    def __repr__(self):
        return '<SamFile Object filename:{0}>'.format(self.filename)


def _parse_line(samline):
    items = samline.split(TAB)
    (qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen,
     seq, qual) = items[:11]
    tags = items[11:]
    return Sam(qname, flag, rname, pos, mapq, cigar, rnext, pnext,
               tlen, seq, qual, *tags)


def parse(samfile):
    with open(samfile, 'r') as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            line = line.rstrip()
            if not line:
                continue
            yield _parse_line(line)


def read(samfile):
    header = []
    handle = open(samfile, 'r')
    offset = 0
    while True:
        line = handle.readline()
        if line.startswith('@'):
            header.append(line.rstrip())
            offset = handle.tell()
        else:
            break

    handle.seek(offset, 0)
    return SamFile(samfile, header, handle, offset)

