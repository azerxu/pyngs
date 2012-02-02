#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: libsff.py
# This file reference from <GSFLX Data Analysis Software Manual.pdf(Page 529)>
# and index part reference from bioPython SffIO.py
# **********************************************************************

"""The Standard Flowgram File is used to store the information on one or many
GS FLX reads and their trace data. Sequencing reads obtained using the Genome
Sequencer FLX System differ from reads obtained using more traditional methods
('Sanger sequencing') in that the GS FLX data does not provide individual base
measurements from which basecalls can be derived. Instead, it provides measure
ments that estimate the length of each homopolymer stretch in the sequence
(e.g. in 'AAATGG', 'AAA' is a 3-mer stretch of A's, 'T' is a 1-mer stretch of
T's and 'GG' is a 2-mer stretch of G’s). A basecalled sequence is then derived
by converting each estimate into a homopolymer stretch of that length and
concatenating the homopolymers.
The .sff file format consists of three sections:
    a common header section occurring once in the file;
    then for each read stored in the file, a read header section and a read
    data section.
    The data in each section consists of a combination of numeric and character
    data; the specific fields for each section are defined below. The sections
    adhere to the following rules:
      I. The standard Unix types uint8_t, uint16_t, uint32_t and uint64_t are
         used to define 1, 2, 4 and 8 byte numeric values.
     II. All multi-byte numeric values are stored using big endian byteorder
         (same as the SCF file format).
    III. All character fields use single-byte ASCII characters.
     IV. Each section definition ends with an 'eight_byte_padding' field, which
         consists of 0 to 7 bytes of padding, so that the byte length of each
         section is divisible by 8 (and hence the next section is aligned on an
         8-byte boundary).
"""

# The standard Unix types uint8_t, uint16_t, uint32_t and uint64_t are
# used to define 1, 2, 4 and 8 byte numeric values.
# char       byte                1        >c
# char[]     string              1        >s
# uint8_t    unsigned char       1        >B
# uint16_t   unsigned short      2        >H
# uint32_t   unsigned int        4        >I
# uint64_t   unsigned long long  8        >Q

__version__ = "0.0.4"
__author__ = "azer.xu"
__email__ = "azer.xu@gmail.com"
__credits__ = ["azer.xu"]
__license__ = "MIT"

import struct
import sys
import getopt

ASCII_A = 65                            # ascii A number
ASCII_0 = 48                            # ascii 0 number

# Format
_SFF = 0x2E736666                       # SFF magic number, as '.sff'
_MFT = 0x2E6D6674                       # MTF Index magic number,as '.mft'


def base36(num, base=36):
    """Calculate number to a base-36 encoding string
    values 0-25 are the letters ‘A’-‘Z’ and
    values 26-35 are the digits ‘0’-‘9’
    Arguments:
    - `num`: a Decimal Number
    """
    _str = []
    while num:
        num, mod = num // base, num % base
        if mod < 26:
            _str.append(chr(ASCII_A + mod))
        else:
            _str.append(chr(ASCII_0 + mod - 26))
    return ''.join(_str[::-1])


def debase36(s):
    """Decode Base36 string to number
    string char 'A'-'Z' are 0 - 25 and
    string char '0'-'9' are 26 - 35
    """
    total = 0
    for c in s:
        num = ord(c)
        if c.isdigit():
            num = num - ASCII_0 + 26
        elif c.isalpha():
            num = num - ASCII_A
        else:
            raise ValueError('Illegal Base36 char: {0}'.format(c))
        total = total * 36 + num
    return total


def timestamp(year=2000, month=1, day=1, hour=1, minute=1, second=1):
    """The timestamp is encoded by computing a “total” value as shown below,
    then converting it into a base-36 string:
        total =
           (year - 2000) * 13 * 32 * 24 * 60 * 60 +
           month * 32 * 24 * 60 * 60 +
           day * 24 * 60 * 60 +
           hour * 60 * 60 +
           minute * 60 +
           second;
    """
    total = (year - 2000) * 13 * 32 * 24 * 60 * 60 + \
            month * 32 * 24 * 60 * 60 + \
            day * 24 * 60 * 60 + \
            hour * 60 * 60 + \
            minute * 60 + \
            second
    return base36(total)


def detimestamp(date):
    """Decode timestamp string to year, month, day, hour, minute, sec
       total =
           (year - 2000) * 13 * 32 * 24 * 60 * 60 +
           month * 32 * 24 * 60 * 60 +
           day * 24 * 60 * 60 +
           hour * 60 * 60 +
           minute * 60 +
           second;
    """
    total = debase36(date)
    total, sec = total // 60, total % 60
    total, minute = total // 60, total % 60
    total, hour = total // 24, total % 24
    total, day = total // 32, total % 32
    year, month = total // 13, total % 13
    year += 2000
    return year, month, day, hour, minute, sec


def random_hash(name):
    """randomizing “hash” character to enhance uniqueness
    Since two Runs may be started at the same second, an additional
    base-36 character is generated by hashing the full rigRunName to a
    base-31 number (the highest prime below 36)
    Arguments:
    - `name`: rig_run_name
              eg. R_yyyy_mm_dd_hh_min_sec_machineName_userName_uniqueRunName
    """
    chval = 0
    for c in name:
        chval += ord(c)
        chval %= 31
    ch = chr(ASCII_A + chval) if chval < 26 else chr(ASCII_0 + chval - 26)
    return ch


def xyloc(x, y):
    """The X,Y location is encoded by computing a total value of
    “X * 4096 + Y” and encoding that as a five character, base-36 string
    Arguments:
    - `x`: x-axis coordinary
    - `y`: y-axis coordinary
    """
    return base36(x * 4096 + y)


def dexyloc(s):
    """Decode xy location string to x, y value
    """
    xy = debase36(s)
    x, y = xy // 4096, xy % 4096
    return x, y


def parse_name(name):
    """This Doc is reference by GSFLX Data Analysis Software Manual
    From page 527-580
    
    454 “Universal” Accession Numbers
    
    The standard 454 read identifiers, used in Genome Sequencer System data
    analysis software versions prior to 1.0.52 (early GS 20), have the format
    “rank_x_y” (as in 003048_1034_0651), where “rank” is a ranking of the well
    in a region by signal intensity, and “x” and “y” are the pixel location of
    the well’s center on the sequencing Run images. This identifier is
    guaranteed to be unique only within the context of a single sequencing
    Run, and may or may not be unique across specific sets of Runs.
    
    To allow for the combination of reads across larger data sets, a more
    unique accession number format has been developed.
    
    An accession in this format is a 14 character string, as in
    C3U5GWL01CBXT2, and consist of 4 components:
        C3U5GW - a six character encoding of the timestamp of the Run
        L      - a randomizing “hash” character to enhance uniqueness
        01     - the region the read came from, as a two-digit number
        CBXT2  - a five character encoding of the X,Y location of the well

    The timestamp, hash character and X,Y location use a base-36 encoding
    (where values 0-25 are the letters ‘A’-‘Z’ and the values 26-35 are the
    digits ‘0’-‘9’). An accession thus consists only of letters and digits,
    and is case-insensitive.
    """
    stamp = name[:6]
    rand_hash = name[6]
    region = name[7:9]
    xy = name[9:]
    date = detimestamp(stamp)           # year, month, day, hour, minute, sec
    region = int(region)
    x, y = dexyloc(xy)
    return date, region, (x, y)


class _CommonHeader(object):
    """The common header section consists of the following fields:
    Field name       Format        Properties
    --------------------------------------------------------------------------
    magic_number     uint32_t      The magic_number field value is 0x2E736666,
                                   the uint32_t encoding of the string '.sff'.
    
    version          char[4]       The version number is 0001, or the byte
                                   array '\0\0\0\1'.
    
    index_offset     uint64_t      The index_offset and index_length fields
    index_length     uint32_t      are the offset and length of an optional
                                   index of the reads in the SFF file. If no
                                   index is included in the file, both fields
                                   must be 0.

    number_of_reads  uint32_t      The number_of_reads field should be set to
                                   the number of reads stored in the file.

    header_length    uint16_t      The header_length field should be the total
                                   number of bytes required by this set of
                                   header fields, and should be equal to '31 +
                                   number_of_flows_per_read + key_length',
                                   rounded up to the next value divisible by 8

    key_length       uint16_t      The key_length field should be set to the
                                   length of the key sequence used for these
                                   reads

    number_of_flows  uint16_t      The number_of_flows_per_read should be set
    _per_read                      to the number of flows for each of the reads
                                   in the file.

    flowgram_format_ uint8_t       The flowgram_format_code should be set to
    code                           the format used to encode code each of the
                                   flowgram values for each read. Currently,
                                   only one flowgram format has been adopted,
                                   so this value should be set to 1. The
                                   flowgram format code 1 stores each value as
                                   a uint16_t, where the floating point
                                   flowgram value is encoded as (int) round(
                                   value * 100.0)', and decoded as
                                   (storedvalue / 100.0). In other words, the
                                   values are stored as an integer encoding of
                                   a limited precision floating point value,
                                   keeping 2 places to the right of the
                                   decimal point, and capping the values at
                                   655.35.

    flow_chars       char[number_  The flow_chars should be set to the array
                     _flows_per_   of nucleotide bases (A, C, G‚ T) that
                     read]         correspond to the nucleotides used for each
                                   flow of each read. The length of the array
                                   should equal number_of_flows_per_read. Note
                                   that the flow_chars field is not
                                   null-terminated.

    key_sequence    char[key_      The key_sequence field should be set to
                    length]        the nucleotide bases of the key sequence
                                   used for these reads. Note that the
                                   key_sequence field is not null-terminated.

    eight_byte_     uint8_t[*]     If any eight_byte_padding bytes exist in
    padding                        the section, they should have a byte value
                                   of 0.
    --------------------------------------------------------------------------
    """
    def __init__(self, magic, version, idx_offset, idx_len, reads_num,
                 hlen, klen, flen, ffmt, flows, keys):
        self.magic = magic              # sff magic number
        self.version = version          # sff version
        self.idx_offset = idx_offset    # sff index offset if exist
        self.idx_len = idx_len          # sff index length if exist
        self.reads_num = reads_num      # sff reads number
        self.hlen = hlen                # sff common header length
        self.klen = klen                # sff key length
        self.flen = flen                # sff flow number for per read
        self.ffmt = ffmt                # flowgram format code
        self.flows = flows              # flowgram flow_chars
        self.keys = keys                # key sequence


class _ReadHeader(object):
    """Read Header Info
    Field name       Format         Properties
    --------------------------------------------------------------------------
    read_header_     uint16_t       The read_header_length should be set to the
    length                          length of the read header for this read,
                                    and should be equal to '16 + name_length'
                                    rounded up to the next value divisible by 8

    name_length      uint16_t       The name_length field should be set to the
                                    length of the read‘s accession or name.

    number_of_bases  uint32_t       The number_of_bases should be set to the
                                    number of bases called for this read.

    clip_qual_left   uint16_t       The clip_qual_left and clip_adapter_left
                                    fields should be set to the position of the
                                    first base after the clipping point, for
                                    quality and/or an adapter sequence, at the
                                    beginning of the read. If only a combined
                                    clipping position is computed, it should be
                                    stored in clip_qual_left.
    clip_qual_       uint16_t
    right                           The clip_qual_right and clip_adapter_right
                                    fields should be set to the position of the
                                    last base before the clipping point, for
                                    quality and/or an adapter sequence, at the
                                    end of the read. If only a combined
                                    clipping position is computed, it should be
                                    stored in clip_qual_right.
    clip_adapter_    uint16_t
    left                            Note that the position values use 1-based
                                    indexing, so the first base is at position
                                    1. If a clipping value is not computed, the
                                    field should be set to 0. Thus, the first
    clip_adapter_    unit16_t       base of the insert is:
    right                                max(1, max(clip_qual_left,
                                             clip_adapter_left))
                                    ...and the last base of the insert is:
                                    min((clip_qual_right == 0 ? number_of_bases
                                          : clip_qual_right),
                                         (clip_adapter_right == 0 ?
                                         number_of_bases : clip_adapter_right))

    name             char[name_     The name field should be set to the string
                     length]        of the read‘s accession or name. Note that
                                    the name field is not null-terminated.

    eight_byte_      uint8_t[*]     If any eight_byte_padding bytes exist in
    padding                         the section, they should have a byte value
                                    of 0.
    --------------------------------------------------------------------------
    """
    def __init__(self, header_length, name_length, base_num, qual_left,
                 qual_right, adapter_left, adapter_right, name):
        self.hlen = header_length
        self.nlen = name_length
        self.bnum = base_num
        self.qleft = qual_left - 1      # make 1-base value to 0-base value
        self.qright = qual_right
        self.aleft = adapter_left - 1   # make 1-base value to 0-base value
        self.aright = adapter_right
        self.name = name


class _ReadData(object):
    """The read data section consists of the following fields:
    Field name       Format         Properties
    --------------------------------------------------------------------------
    flowgram_values  uint*_t[number The flowgram_values field contains the
                    _of_flows]      homopolymer stretch estimatesfor each flow
                                    of the read. The number of bytes used for
                                    each value depends on the common header
                                    flowgram_format_code value (where the
                                    current value uses a uint16_t for each
                                    value).
    
    flow_index_per   uint8_t[number The flow_index_per_base field contains the
    _base            _of_bases]     flow positions for each base in the called
                                    sequence (i.e., for each base, the position
                                    in the flowgram whose estimate resulted in
                                    that base being called). Note that these
                                    values are “incremental” values, i.e. the
                                    stored position is the offset from the
                                    previous flow index in the field. All
                                    position values (prior to their incremental
                                    encoding) use 1-based indexing, so the
                                    first flow is flow 1.
    
    Bases            char[number_   The bases field contains the basecalled
                     of_ bases]     nucleotide sequence.
    
    quality_scores   uint8_t[number The quality_scores field contain the
                    _of_bases]      quality scores for each of the bases in
                                    the sequence, where the values use the
                                    standard -log10 probability scale.
    
    eight_byte_     uint8_t[*]      If any eight_byte_padding bytes exist in
    padding                         the section, they should have a byte value
                                    of 0.
    --------------------------------------------------------------------------
    """
    __slots__ = ['bases', 'quals', 'flows', 'idxs']
    
    def __init__(self, flows, idxs, bases, quals):
        self.bases = bases              # basecalled nucleotide sequence
        self.quals = quals              # quality scores for each base in seq
        self.flows = flows              # flowgram_values estimates each flow
        self.idxs = idxs                # base index of flowgram


class SffRead(object):
    """Recoder Sff Read info
    """
    
    def __init__(self, header, data):
        self.header = header            # _ReadHeader Instance
        self.data = data                # _ReadData Instance

    def __len__(self):
        """return the length of all bases called by Machine
        """
        return self.header.bnum

    @property
    def lclip(self):
        """return left clip
        base of the first insert (left) is:
            max(1, max(clip_qual_left, clip_adapter_left))
        """
        # use 0 because of python is 0-base not 1-base
        return max((0, self.header.qleft, self.header.aleft))

    @property
    def rclip(self):
        """return right clip
        the last base of the insert (right) is:
        min( (clip_qual_right == 0 ? number_of_bases : clip_qual_right),
             (clip_adapter_right == 0 ? number_of_bases : clip_adapter_right)
           )
        """
        return min((self.header.qright if self.header.qright else \
                    self.header.bnum,
                    self.header.aright if self.header.aright else \
                    self.header.bnum))

    @property
    def bases(self):
        """return called all bases
        """
        return self.data.bases

    @property
    def flows(self):
        """return all flowgrams read by Machine
        """
        return ' '.join(('{0:.2f}'.format(_flow) for _flow in self.data.flows))

    @property
    def idxs(self):
        """return index of each base called in flows
        """
        return ' '.join(('{0}'.format(_idx) for _idx in self.data.idxs))

    @property
    def quals(self):
        """return reads quality for each bases
        """
        return ' '.join(('{0}'.format(_qual) for _qual in self.data.quals))

    @property
    def length(self):
        """Reads bases length after trimed
        """
        return self.rclip - self.lclip

    @property
    def full_length(self):
        """Reads bases length with all bases
        """
        return self.header.bnum

    @property
    def date(self):
        """return (year, month, day, hour, minute, sec)
        to indicate when the reads generate
        """
        return detimestamp(self.name[:6])

    @property
    def region(self):
        """return read in which region
        """
        return int(self.name[7:9])

    @property
    def xyloc(self):
        """return read coordinary (x, y) in plate
        """
        return dexyloc(self.name[9:])

    @property
    def name(self):
        """Return reads name
        """
        return self.header.name
    
    @property
    def title(self):
        """Generate a title to contain reads info
        """
        x, y = self.xyloc
        return '{0} x={1} y={2} region={3} length={4}'.format(
            self.name, x, y, self.region, self.length)

    @property
    def full_title(self):
        """Generate a title to contain reads info
        """
        # date = "{0}_{1}_{2}_{3}_{4}_{5}".format(*self.date)
        x, y = self.xyloc
        return '{0} x={1} y={2} region={3} length={4} lclip={5} rclip={6}'\
               .format(self.name, x, y, self.region, self.full_length,
                       self.lclip + 1, self.rclip)

    @property
    def full_bases(self):
        """return trimed seq in lower case, untrimed in upper case"""
        return ''.join((self.bases[:self.lclip].lower(),
                         self.bases[self.lclip:self.rclip],
                         self.bases[self.rclip:].lower()))

    @property
    def fasta(self):
        """return a fasta format record
        """
        return '\n'.join(
            ('>{0}'.format(self.title), self.bases[self.lclip:self.rclip]))

    @property
    def full_fasta(self):
        """return a full fasta format record
        """
        return '\n'.join(('>{0}'.format(self.full_title), self.full_bases))

    @property
    def phred33_quals(self):
        return ''.join((chr(33 + qual) for qual in self.data.quals))

    @property
    def phred64_quals(self):
        return ''.join((chr(64 + qual) for qual in self.data.quals))

    @property
    def fastq(self):
        """return a cliped fastq format record
        """
        return '\n'.join(('@{0}'.format(self.title),
                          self.bases[self.lclip:self.rclip],
                          '+', self.phred33_quals[self.lclip:self.rclip]))

    @property
    def full_fastq(self):
        """return a full fastq format record
        """
        return '\n'.join(
            ('@{0}'.format(self.full_title), self.full_bases, '+',
             self.phred33_quals))

    @property
    def flow(self):
        """return name/flows as a tab-delimited line
        """
        return '\t'.join((self.name, self.flows))

    @property
    def qual(self):
        """return fasta format quality record
        """
        return '\n'.join(('>{0}'.format(self.title),
                          self.quals[self.lclip:self.rclip]))

    @property
    def full_qual(self):
        """
        """
        return '\n'.join(('>{0}'.format(self.full_title), self.quals))

    @property
    def seq(self):
        """return name/seq as tab-delimited line
        """
        return '\t'.join((self.name, self.bases[self.lclip:self.rclip]))

    @property
    def full_seq(self):
        return '\t'.join((self.name, self.bases))

    @property
    def tab(self):
        """return the seq/qual/flow as tab-delimited line
        """
        return '\t'.join(self.name, self.bases, self.quals,
                         self.flows)

    def __repr__(self):
        return '\n'.join((self.title, self.flows, self.idxs,
                          self.bases, self.quals))


class _Index(object):
    """Recorder any existing Roche style XML meta data
    """
    def __init__(self, xml=''):
        self.xml = xml


class SffFile(object):
    """Record Sff file info
    """
    def __init__(self, fname):
        self._handle = open(fname, 'rb')
        self.header = self.__parse_common_header()
        self.index = self.__parse_index()
        self._handle.seek(self.header.hlen)

    def _read(self, num):
        return self._handle.read(num)

    def _get(self, fmt):
        fmt_size = struct.calcsize(fmt)
        return struct.unpack(fmt, self._read(fmt_size))

    def __parse_common_header(self):
        self._handle.seek(0)
        fmt = '>I4BQIIHHHB'
        assert struct.calcsize(fmt) == 31
        (magic, v1, v2, v3, v4, index_offset, index_length,
         number_of_reads, header_length, key_length,
         number_of_flows_per_read,
         flowgram_format_code) = struct.unpack(fmt, self._read(31))
        assert magic == _SFF, "Not Legal sff file"
        
        version = '{0}{1}{2}{3}'.format(v1, v2, v3, v4)
        
        flow_chars = struct.unpack('>{0}s'.format(number_of_flows_per_read),
                                   self._read(number_of_flows_per_read))
        
        key_sequence = struct.unpack('>{0}s'.format(key_length),
                                     self._read(key_length))
        
        return _CommonHeader(
            magic, version, index_offset, index_length, number_of_reads,
            header_length, key_length, number_of_flows_per_read,
            flowgram_format_code, flow_chars, key_sequence)

    def __parse_index(self):
        """Parse Roche 454 Index xml meta data and read index
        """
        idx_len = self.header.idx_len
        idx_offset = self.header.idx_offset
        if not idx_offset or not idx_len:
            return Index()
        self._handle.seek(idx_offset)
        fmt1 = '>I4s'                   # magic and version format
        fmt1_size = struct.calcsize(fmt1)
        
        (idx_magic, idx_version) = struct.unpack(fmt1, self._read(fmt1_size))
        
        assert idx_magic == _MFT, "Illegal Roche 454 Index magic Number"
        
        fmt2 = '>LL'                    # xml size and data size
        fmt2_size = struct.calcsize(fmt2)
        
        (xml_size, data_size) = struct.unpack(fmt2, self._read(fmt2_size))
        
        xml = self._handle.read(xml_size)
        
        # **************************************************
        # # index data, but don't know how to use
        # fmt3 = '>14s2xI'
        # fmt3_size = struct.calcsize(fmt3)
        # for i in xrange(self.header.reads_num):
        #     (acc, loc) = struct.unpack(fmt3, self._handle.read(fmt3_size))
        #     print acc, loc
        #     break
        # **************************************************
        return _Index(xml)

    def __parse_read(self):
        """parse read header and read data, return a SffRead instance
        """
        # ********** Start Read header Block **********
        # H read_header_length
        # H name_length
        # I number_of_bases
        # H clip_qual_left
        # H clip_qual_right
        # H clip_adapter_left
        # H clip_adapter_right
        fmt = ">HHIHHHH"                # reads header format
        fmt_size = struct.calcsize(fmt)
        # assert fmt_size == 16, "Reads Header length is not equal 16"
        
        # get read header info
        (read_header_length, name_length, base_num, clip_qual_left,
         clip_qual_right, clip_adapter_left,
         clip_adapter_right) = struct.unpack(fmt, self._read(fmt_size))
        
        name = self._read(name_length)   # get read name
        
        # calc heading padding make data aligned with 8
        padding = read_header_length - 16 - name_length
        if padding:
            self._read(padding)
        # ********** Finish Read Header Block **********
        
        _read_header = _ReadHeader(     # generate ReadHeader instance
            read_header_length, name_length, base_num, clip_qual_left,
            clip_qual_right, clip_adapter_left, clip_adapter_right, name)
        
        # ********** Start Read Data Block **********
        flow_fmt = '>{0}H'.format(self.header.flen)
        flows = self._get(flow_fmt)
        flows = [flow / 100.0 for flow in flows]
        
        idx_fmt = '>{0}B'.format(base_num)
        idxs = self._get(idx_fmt)
        
        base_fmt = '>{0}s'.format(base_num)
        (bases,) = self._get(base_fmt)
        
        qual_fmt = '>{0}B'.format(base_num)
        quals = self._get(qual_fmt)
        
        padding = (2 * self.header.flen + 3 * base_num) % 8
        if padding:
            self._read(8 - padding)
        # ********** Finish Read Data Block **********
        
        _read_data = _ReadData(flows, idxs, bases, quals)
        
        return SffRead(_read_header, _read_data)

    @property
    def mft(self):
        """return index xml info"""
        return self.index.xml

    @property
    def reads(self):
        """return each read one by one from sff file"""
        for i in xrange(self.header.reads_num):
            yield self.__parse_read()

    @property
    def readnum(self):
        """return how many reads in sff file"""
        return self.header.reads_num


def read(fname):
    """Read Sff file in rb mode and return a Iterator of SffReads
    """
    sff = SffFile(fname)
    return sff.reads


def errlog(msg):
    print >>sys.stderr, msg


def show_verson():
    print "Mysff version: {0}".format(__version__)
    exit()


def show_usage():
    print """Mysff is Roche 454 Sff file extractor (Version: {0})
    
    Usage: mysff [opts] sfffile1 sfffile2 sfffile3 ...
    Opts:
       -a or --accno    Output just the accessions
       -s or --seq      Output just the sequences
       -q or --qual     Output just the quality scores
       -f or --flow     Output just the flowgrams
       -t or --tab      Output the seq/qual/flow as tab-delimited lines
       -n or --notrim   Output the untrimmed sequence or quality scores
       -m or --mft      Output the manifest text
       -r or --readnum  Output the total reads number
       -u or --fastq    Output the fastq format (quality is phred33 as default)
       -v or --version  Output the mysff current version
       -h or --help     Output help (show this message)

    if no opts given, then fastq is set, trim is set to generate trimed fastq.

    sfffile is full path file name, if sff file name don't contain Full path,
    then the file is looked in current folder. mysff can extract multi sff
    file. the opt -n or --notrim only affect -s -q -u (--seq, --qual, --fastq).
    if -n or --notrim is set, the sequence between left trim and right trim in
    upper case, and the sequence before left trim and after right trim in lower
    case.
    """.format(__version__)
    exit()                              # quit program


def main(*argv):
    if not len(argv):                   # no args given
        show_usage()

    try:
        optlst, args = getopt.getopt(
            argv, 'asqftnmuvhr', ['accno', 'seq', 'qual', 'flow', 'tab',
                                 'notrim', 'mft', 'fastq', 'version', 'help',
                                  'readnum'])

        method = 'fastq'                # set fastq as default extract method
        istrim = True                   # identify trim data or not
        for opt, val in optlst:
            if opt in ('-h', '--help'): # show usage
                show_usage()
            elif opt in ('-v', '--version'): # show version
                show_verson()
            elif opt in ('-t', '--tab'): # output seq/qual/flow
                method = 'tab'
            elif opt in('-m', '--mft'): # output manifest text
                method = 'mft'
            # output untrimmed sequence, quality or fastq
            elif opt in('-n', '--notrim'):
                istrim = False
            elif opt in('-a', '--accno'): # only output accessions
                method = 'name'
            elif opt in ('-s', '--seq'): # output fasta
                method = 'fasta'
            elif opt in ('-f', '--flow'): # output flowgrams
                method = 'flow'
            elif opt in ('-q', '--qual'): # output quality scores
                method = 'qual'
            elif opt in ('-u', '--fastq'): # output fastq
                method = 'fastq'
            elif opt in ('-r', '--readnum'): # output the total reads number
                method = 'readnum'

    except getopt.GetoptError, e:
        errlog("Err: {0}".format(e))
        show_usage()

    # no sff file given
    if not len(args):
        errlog('Err: No Sff file given\n')
        show_usage()

    if not istrim and method in ('fasta', 'fastq', 'qual'):
        method = 'full_' + method

    for arg in args:
        sff = SffFile(arg)
        if method == 'mft':
            print sff.mft
        elif method == 'readnum':
            print arg, sff.readnum
        else:
            for sread in sff.reads:
                print getattr(sread, method)

if __name__ == '__main__':
    main(*sys.argv[1:])
