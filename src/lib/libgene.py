#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: libgene.py
#
# This module is parse Gene structure for annotation SNP
# **********************************************************************

class Exon(object):
    def __init__(self, mark, name, length, strand, start, end):
        self.mark = mark
        self.name = name
        self.length = self.deal_num(length, fix=False)
        self.strand = strand
        self.start = self.deal_num(start, fix=True)
        self.end = self.deal_num(end, fix=False)

    def deal_num(self, x, fix=False):
        # fix means trans x to 0-base or not
        if isinstance(x, str):
            if fix:                     # fix 1-base number to 0-base number
                return int(x) - 1
            else:
                return int(x)
        elif isinstance(x, int):        # if number is int then think
            return x                    # x as 0-base number, no need to fix
        else:
            raise ValueError('Illegal data type: {0}'.format(x))

    def __contains__(self, x):           # used as 'in'
        if self.start <= x < self.end:
            return True
        else:
            return False

    def __repr__(self):
        return '\t'.join((self.mark, self.name, str(self.length),
                          self.strand, str(self.start+1), str(self.end)))


class Gene(Exon):
    def __init__(self, mark, name, length, strand, start, end):
        super(Gene, self).__init__(mark, name, length, strand, start, end)
        self.utr5 = []
        self._exons = []
        self.utr3 = []

    def add(self, exon):
        if exon.name == 'UTR3':
            self.utr3.append(exon)
        elif exon.name == 'UTR5':
            self.utr5.append(exon)
        # elif exon.name.startswith('EXON'):
        #     self._exons.append(exon)
        # elif exon.name in ('CDS', 'rRNA', 'tRNA', 'ncRNA'):
        #     self._exons.append(exon)
        # else:
        #     raise ValueError('Unkown Exon: {0}'.format(repr(exon)))
        else:
            self._exons.append(exon)

    @property
    def exons(self):
        if not self._exons and not self.utr5 and not self.utr3:
            return [Exon('EXON', 'EXON', self.length, self.strand,
                         self.start, self.end)]
        else:
            return self._exons

    @property
    def fullexons(self):
        if self.strand == '+':
            return self.utr5 + self.exons + self.utr3
        elif self.strand == '-':
            return self.utr3 + self.exons + self.utr5

    @property
    def cdslength(self):
        return sum([exon.length for exon in self.exons])

    @property
    def genelength(self):
        return self.length

    @property
    def exoncount(self):
        return len(self.exons)

    @property
    def mrnalength(self):
        return sum([exon.length for exon in self.fullexons])

    def __repr__(self):
        gene = super(Gene, self).__repr__()
        exons = [repr(exon) for exon in self.fullexons]
        return '\n'.join([gene] + exons)


# **********************************************************************
# parse gene list table, the table format as below
# Mark  Name   Length strand  Start   End
#
# An Example:
# GENE  NM109  21124  -  245888  267011
# EXON  EXON7  887    -  245888  246774
# EXON  EXON6  124    -  256157  256280
# EXON  EXON5  225    -  257490  257714
# EXON  EXON4  804    -  258567  259370
# EXON  EXON3  273    -  259805  260077
# EXON  EXON2  16     -  260547  260562
# EXON  EXON1  206    -  266806  267011
# **********************************************************************
def parse(fname, sep='\t'):
    """load gene output table"""
    gene = None

    with open(fname, 'r') as handle:
        for line in handle:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            items = line.split(sep)
            if line.startswith('GENE'):
                if gene:
                    yield gene
                gene = Gene(*items)
            elif line.startswith('EXON'):
                if not gene:
                    raise ValueError(
                        'file structure is Wrong: {0}'.format(line))
                exon = Exon(*items)
                gene.add(exon)
            else:
                raise VauleError('Unsupport Record: {0}'.format(line))

        if gene:
            yield gene


# **********************************************************************
# Below block is parse ucsc record
# **********************************************************************
def _parse_list(lst):
    lst = lst.replace('"', '')
    if lst.endswith(','):
        lst = lst[:-1]
    lst = lst.split(',')
    return [int(each) for each in lst]


def _get_cds_start_split_index(cdsstart, exons):
    for i in xrange(len(exons)):
        if cdsstart < exons[i][0]:
            return i - 1
    return i                        # self.exoncount - 1


def _get_cds_end_split_index(cdsend, exons):
    for i in xrange(len(exons)-1, -1, -1):
        if cdsend > exons[i][-1]:
            return i + 1
    return i                        # 0


def _calc_utr_head(cdsstart, exons):
    if not exons:                       # exons is empty
        return [], []

    utr_idx = _get_cds_start_split_index(cdsstart, exons)

    utr = exons[:utr_idx]

    start, end = exons[utr_idx]         # split exon's start and end

    if cdsstart == start:               # cdsstart equal exon start
        return utr, exons[utr_idx:]

    # split the utr_idx exon to utr and exon
    utr.append((start, cdsstart))

    if cdsstart < end:
        exons = [(cdsstart, end)] + exons[utr_idx+1:]
    else:
        exons = exons[utr_idx+1:]
    return utr, exons


def _calc_utr_tail(cdsend, exons):
    if not exons:                       # exons is empty
        return [], []

    utr_idx = _get_cds_end_split_index(cdsend, exons)

    utr = exons[utr_idx+1:]

    start, end = exons[utr_idx]

    if cdsend == end:
        return utr, exons[:utr_idx+1]

    utr = [(cdsend, end)] + utr

    if cdsend > start:
        exons = exons[:utr_idx] + [(start, cdsend)]
    else:
        exons = exons[:utr_idx]
    return utr, exons

def _calc_exons(cdsstart, cdsend, exonstarts, exonends, strand='+'):
    cdsstart = int(cdsstart)
    cdsend = int(cdsend)
    starts = _parse_list(exonstarts)
    ends = _parse_list(exonends)
    exons = zip(starts, ends)
    utr_head, exons = _calc_utr_head(cdsstart, exons)
    utr_tail, exons = _calc_utr_tail(cdsend, exons)

    mark = 'EXON'
    # parse UTR at head
    name = 'UTR5' if strand == '+' else 'UTR3'
    for start, end in utr_head:
        yield Exon(mark=mark, name=name, length=end-start, strand=strand,
                   start=start, end=end)

    # parse EXONs
    exon_count = len(exons)
    for i, (start, end) in enumerate(exons):
        name = 'EXON{0}'.format(i+1 if strand == '+' else exon_count-i)
        yield Exon(mark=mark, name=name, length=end-start, strand=strand,
                   start=start, end=end)

    # parse UTR at tail
    name = 'UTR3' if strand == '+' else 'UTR5'
    for start, end in utr_tail:
        yield Exon(mark=mark, name=name, length=end-start, strand=strand,
                   start=start, end=end)


def deal_ucsc_record_line(line, sep='\t'):
    try:
        (name, chrom, strand, txstart, txend, cdsstart, cdsend,
         exoncount, exonstarts, exonends, other) = line.split(sep, 10)
    except ValueError, ex:
        print line
        raise ValueError(ex)
    chrom_name = '{0}.{1}'.format(chrom, name)
    txstart = int(txstart)
    txend = int(txend)
    length = txend - txstart
    gene = Gene(mark='GENE', name=chrom_name, length=length, strand=strand,
                start=txstart, end=txend)
    for exon in _calc_exons(cdsstart, cdsend, exonstarts, exonends, strand):
        gene.add(exon)
    return gene


# **********************************************************************
# parse ucsc record and return a Gene instance
# **********************************************************************
def parse_ucsc(fname):
    with open(fname, 'r') as handle:
        for line in handle:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            yield deal_ucsc_record_line(line)


# **********************************************************************
# divide genes to blocks, each block size is 0.1M
# **********************************************************************
def parse_gene_list(fname, block=100000):
    genes = {'-': {}, '+': {}}

    for gene in parse(fname):
        chrom, gname = gene.name.split('.')[-2:]
        strand = gene.strand
        if chrom not in genes[strand]:
            genes[strand][chrom] = {}
        # calc block ids
        start_block_id = gene.start // block
        end_block_id = gene.end // block
        for i in xrange(start_block_id, end_block_id+1):
            if i not in genes[strand][chrom]:
                genes[strand][chrom][i] = []
            genes[strand][chrom][i].append(gene)
    return genes


# **********************************************************************
# test parse module
# **********************************************************************
def test_parse():
    for gene in parse('rat.gene.list'):
        print gene
        print '*' * 70
        print


def test_parse_ucsc():
    for gene in parse_ucsc('rat.refGene.txt'):
        print gene
        # print '*' * 70
        # raw_input()


if __name__ == '__main__':
    import sys
    test_parse_ucsc()

