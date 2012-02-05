#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: test_libmp.py
# **********************************************************************

from pyngs.lib.libmp import run, SENTINEL


def producer(n):
    for i in xrange(n):
        yield i


def consumer(i):
    # print 'pid:', os.getpid(), 'now process', i
    return i, i ** 3


def _reporter(i, j):
    print 'end:', i, j

def reporter(iqueue, nconsumer=1, sentinel=SENTINEL):
    while True:
        item = iqueue.get()
        if item is sentinel:
            nconsumer -= 1
            if nconsumer < 1:           # done
                break
            continue

        _reporter(*item)


if __name__ == '__main__':
    for num in (10000,):
        run(producer=producer, producer_args=(num,),
            consumer=consumer,
            reporter=reporter,
            nconsumer=10
            )

