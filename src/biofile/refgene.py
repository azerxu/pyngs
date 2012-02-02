#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: librefgene.py
#
# This module is parse UCSC refGene.txt file to Gene structure for
# annotation SNP
# **********************************************************************

from libgene import Exon, Gene

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


# CREATE TABLE `refGene` (
#   `bin` smallint(5) unsigned NOT NULL,
#   `name` varchar(255) NOT NULL,
#   `chrom` varchar(255) NOT NULL,
#   `strand` char(1) NOT NULL,
#   `txStart` int(10) unsigned NOT NULL,
#   `txEnd` int(10) unsigned NOT NULL,
#   `cdsStart` int(10) unsigned NOT NULL,
#   `cdsEnd` int(10) unsigned NOT NULL,
#   `exonCount` int(10) unsigned NOT NULL,
#   `exonStarts` longblob NOT NULL,
#   `exonEnds` longblob NOT NULL,
#   `score` int(11) default NULL,
#   `name2` varchar(255) NOT NULL,
#   `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
#   `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
#   `exonFrames` longblob NOT NULL,
def deal_refgene_record_line(line, sep='\t'):
    try:
        (_bin, name, chrom, strand, txstart, txend, cdsstart, cdsend,
         exoncount, exonstarts, exonends, score, name2,
         cdsstartstat, cdsendstat, exonframes) = line.split(sep)
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
def parse(fname):
    with open(fname, 'r') as handle:
        for line in handle:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            yield deal_refgene_record_line(line)


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


def main(args):
    for arg in args:
        for i, gene in enumerate(parse(arg)):
            print gene


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

