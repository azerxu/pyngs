#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: filsam.py
# **********************************************************************
import os
import sys
from pyngs.biofile.sam import read


def get_basename(fname):
    name, ext = os.path.splitext(os.path.basename(fname))
    return name, ext


def get_pair(samfileobj):
    while True:
        try:
            read1 = samfileobj.next()
            read2 = samfileobj.next()
            yield read1, read2
        except StopIteration:
            break


def output_head(samfileobj, out):
    for line in samfileobj.header:
        print >>out, line.rstrip()


def filpair(samfile, nm=2):
    name, ext = get_basename(samfile)
    # reads mapping and unique mapping and mismath less than nm para given
    good = open('{0}.good{1}'.format(name, ext), 'w')

    # both two of the reads not mapping
    unmap = open('{0}.unmap{1}'.format(name, ext), 'w')

    # reads contain soft clip
    soft = open('{0}.soft{1}'.format(name, ext), 'w')

    # one is mapping and the other is not mapping
    onemap = open('{0}.onemap{1}'.format(name, ext), 'w')

    # reads mapping but has more than one location mapping
    repeat = open('{0}.rep{1}'.format(name, ext), 'w')

    # one of reads mapped in N region
    nregion = open('{0}.n{1}'.format(name, ext), 'w')

    # one of reads mapped as mate sw method
    matesw = open('{0}.mate{1}'.format(name, ext), 'w')

    # reads has uniq mapping but has multi suboptional alignment
    multi = open('{0}.mul{1}'.format(name, ext), 'w')

    # reads unique mapping but has more error in reads
    miserr = open('{0}.nm{1}{2}'.format(name, nm, ext), 'w')

    # reads mapping to two chromas or scaffold
    cross = open('{0}.cross{1}'.format(name, ext), 'w')

    # reads filter error
    err = open('{0}.err{1}'.format(name, ext), 'w')

    samfileobj = read(samfile)

    output_head(samfileobj, good)

    for read1, read2 in get_pair(samfileobj):
        try:
            if read1.is_unmapped and read2.is_unmapped:
                # both two reads not mapped
                print >>unmap, read1
                print >>unmap, read2
            elif read1.is_unmapped or read2.is_unmapped:
                # one is mapped, other not mapped
                print >>onemap, read1
                print >>onemap, read2
            elif not read1.mapq or not read2.mapq: # mapped in repeat region
                print >>repeat, read1
                print >>repeat, read2
            elif read1.xt == 'N' or read2.xt == 'N': # map in N block
                print >>nregion, read1
                print >>nregion, read2
            elif read1.rnext != read2.rnext: # reads mapping to two scaffold
                print >>cross, read1
                print >>cross, read2
            elif 'S' in read1.cigar or 'S' in read2.cigar: # soft clip in reads
                print >>soft, read1
                print >>soft, read2
            elif read1.xt == 'M' or read2.xt == 'M': # one read is use bwasw
                print >>matesw, read1
                print >>matesw, read2
            elif (read1.x0 > 1 or read2.x0 > 1 or
                  read1.x0+read1.x1 > 1 or read2.x0+read2.x1 > 1):
                # has multi mapping position
                print >>multi, read1
                print >>multi, read2
            elif read1.nm > nm or read2.nm > nm: # too much mismatchs
                print >>miserr, read1
                print >>miserr, read2
            else:                       # both reads fit the requirement
                print >>good, read1
                print >>good, read2
        except:
            # unexpected results
            print >>err, read1
            print >>err, read2

    for out in (good, unmap, repeat, multi, matesw, miserr,
                err, onemap, cross, nregion):
        out.close()


if __name__ == '__main__':
    map(filpair, sys.argv[1:])

