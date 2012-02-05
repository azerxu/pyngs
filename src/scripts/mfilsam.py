#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: mfilsam.py
#
# use libmp.py module to filter sam
# **********************************************************************
import os
import sys
from pyngs.biofile.sam import Sam
from pyngs.lib.libmp import run, SENTINEL
import getopt


# define CONSTANCE
GOOD = 0
UNMAP = 1
ONEMAP = 2
REPEAT = 3
CROSS = 4
SOFT = 5
NREGION = 6
MULTI = 7
MATESW = 8
MISERR = 9
ERR = 10


def parse_line(line):
    items = line.split('\t')
    return Sam(*items)


def get_basename(fname):
    name, ext = os.path.splitext(os.path.basename(fname))
    return name, ext


def get_pair(samfile):                  # producer
    with open(samfile, 'r') as handle:
        while True:
            line1 = handle.readline().rstrip()
            if not line1.startswith('@'):
                break

        line2 = handle.readline().rstrip()
        if line1 and line2:
            yield line1, line2
        else:
            return

        while True:
            line1 = handle.readline().rstrip()
            line2 = handle.readline().rstrip()
            if line1 and line2:
                yield line1, line2
            else:
                break


def classify_pair(line1, line2, nm=2):  # consumer
    read1 = parse_line(line1)
    read2 = parse_line(line2)
    try:
        if read1.is_unmapped and read2.is_unmapped:
            # both two reads not mapped
            mark = UNMAP
        elif read1.is_unmapped or read2.is_unmapped:
            # one is mapped, other not mapped
            mark = ONEMAP
        elif read1.xt == 'N' or read2.xt == 'N': # map in N block
            mark = NREGION
        elif not read1.mapq or not read2.mapq: # mapped in repeat region
            mark =  REPEAT
        elif read1.rname != read2.rname:
            # reads mapping to two scaffold
            mark = CROSS
        elif 'S' in read1.cigar or 'S' in read2.cigar: # soft clip in reads
            mark = SOFT
        elif read1.xt == 'M' or read2.xt == 'M': # one read is use bwasw
            mark = MATESW
        elif (read1.x0 > 1 or read2.x0 > 1 or
              read1.x0+read1.x1 > 1 or read2.x0+read2.x1 > 1):
            # has multi mapping position
            mark = MULTI
        elif read1.nm > nm or read2.nm > nm: # too much mismatchs
            mark = MISERR
        else:                       # both reads fit the requirement
            mark = GOOD
    except:
        mark = ERR                      # unexpected results
    return mark, line1, line2


# reporter
def output_pair(iqueue, samfile, nm=2, nconsumer=1, sentinel=SENTINEL):
    name, ext = get_basename(samfile)
    # reads mapping and unique mapping and mismath less than nm para given
    good = open('{0}.good{1}'.format(name, ext), 'w')

    # both two of the reads not mapping
    unmap = open('{0}.unmap{1}'.format(name, ext), 'w')

    # one is mapping and the other is not mapping
    onemap = open('{0}.onemap{1}'.format(name, ext), 'w')

    # reads mapping but has more than one location mapping
    repeat = open('{0}.rep{1}'.format(name, ext), 'w')

    # one of reads mapped in N region
    nregion = open('{0}.n{1}'.format(name, ext), 'w')

    # reads contain soft clip
    soft = open('{0}.soft{1}'.format(name, ext), 'w')

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

    outs = (good, unmap, onemap, repeat, cross, soft, nregion, multi, matesw,
            miserr, err)

    while True:
        item = iqueue.get()
        if item is sentinel:            # check item is sentinel or not
            nconsumer -= 1              # if nconsumer all done their jobs
            if nconsumer < 1:           # then break the circle
                break
            continue
        mark, read1, read2 = item
        out = outs[mark]
        print >>out, read1
        print >>out, read2

    # all done
    for out in outs:
        out.close()


def show_usage():
    print 'Usage: mfilsam.py [-nm] samfil1 samfile2 ...'
    print '       nm default value is 2'
    exit()


def main(argv):
    try:
        optlst, args = getopt.getopt(
            argv, 'ht:n:c:', ['help', 'pnum', 'nm', 'nconsumer'])
        nm = 2
        pnum = 2
        nconsumer = 2
        for opt, val in optlst:
            if opt in ('-h', '--help'): # show help message
                show_usage()
            elif opt in ('-t', '--pnum'): # process number
                pnum = int(val)
            elif opt in ('-n', '--nm'): # mismatch num max value
                nm = int(val)
            elif opt in ('-c', '--nconsumer'): # consumer process number
                nconsumer = int(val)
            else:
                show_usage()
    except GetoptError:
        show_usage()

    if not args:                        # no samfile given
        show_usage()

    for arg in args:                    # run each file separate
        run(producer=get_pair, producer_args=(arg,),
            consumer=classify_pair, consumer_kwargs=dict(nm=nm),
            reporter=output_pair, reporter_args=(arg,),
            reporter_kwargs=dict(nm=nm),
            sentinel=SENTINEL, nconsumer=nconsumer, pnum=pnum)


if __name__ == '__main__':
    main(sys.argv[1:])

