#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: test_libmp.py
# **********************************************************************

from pyngs.lib.libmp import run
# import time
# import os


def producer(n):
    for i in xrange(n):
        yield i


def consumer(i):
    # print 'pid:', os.getpid(), 'now process', i
    return i, i ** 3


def reporter(i, j):
    print 'end:', i, j


if __name__ == '__main__':
    for num in (10000,):
        run(producer=producer, producer_args=(num,),
            consumer=consumer,
            reporter=reporter,
            nconsumer=10
            )

