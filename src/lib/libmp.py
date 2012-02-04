#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: libmp.py
# **********************************************************************
"""libmp module for run such mode:
        | consumer |
producer| consumer | reporter
        | consumer |
        |   ....   |
"""


from multiprocessing import Process, Queue

DONE = None


def _get_producer(target=None, name='producer', args=(), kwargs={},
                  oqueue=None, nconsumer=1):
    def _func(*args, **kwargs):
        for product in target(*args, **kwargs):
            oqueue.put(product)

        for cid in range(nconsumer):
            oqueue.put(DONE)

    return Process(target=_func, name=name, args=args, kwargs=kwargs)


def _get_consumer(target=None, name='consumer', args=(), kwargs={},
                  iqueue=None, oqueue=None, nproducer=1):

    def _func(*args, **kwargs):
        _counter = nproducer
        while True:
            item = iqueue.get()
            if item is DONE:
                _counter -= 1
                if _counter < 1:        # check producer is done or not
                    break               # if yes, then break the circle
                continue                # if not wait to producer finished

            if isinstance(item, tuple) or isinstance(item, list):
                _args = item + args
            else:
                _args = (item,) + args
            res = target(*_args, **kwargs)
            oqueue.put(res)

        oqueue.put(DONE)

    return Process(target=_func, name=name, args=args, kwargs=kwargs)


def _get_reporter(target=None, name='reporter', args={}, kwargs={},
                  iqueue=None, nconsumer=1):

    def _func(*args, **kwargs):
        _counter = nconsumer
        while True:
            item = iqueue.get()
            if item is None:
                _counter -= 1
                if _counter < 1:        # check consumer is finished or not
                    break               # break the circle
                continue                # if not finished then wait to finish

            if isinstance(item, tuple) or isinstance(item, list):
                _args = item + args
            else:
                _args = (item,) + args
            target(*_args, **kwargs)

    return Process(target=_func, name=name, args=args, kwargs=kwargs)


def run(producer=None, producer_name='producer', producer_args=(),
        producer_kwargs={}, consumer=None, consumer_name='consumer',
        consumer_args=(), consumer_kwargs={}, reporter=None,
        reporter_name='reporter', reporter_args=(), reporter_kwargs={},
        pnum=None, nconsumer=None):
    # first check producer, consumer and reporter method setted
    if not producer or not consumer or not reporter:
        raise ValueError("producer, consumer, reporter can't be empty")

    # construct pipe to transport data
    producer_queue = Queue()
    consumer_queue = Queue()

    # producer number always set to 1 in this version
    nproducer = 1

    # reporter number always set to 1 in this version
    nreporter = 1

    # calc nconsumer if not given
    if not nconsumer:
        if not pnum:                    # if both nconsumer and process set
            nconsumer = 1               # then set consumer number to 1
        else:
            # set nconsumer as total process minus nproducer and nreporter
            nconsumer = pnum - nproducer - nconsumer

    if nconsumer < 1:                   # if nconsumer less than 1
        nconsumer = 1                   # then set it to default 1

    # create producer Process
    _producer = _get_producer(target=producer,name=producer_name,
                              args=producer_args, kwargs=producer_kwargs,
                              nconsumer=nconsumer, oqueue=producer_queue)
    _producer.start()                   # start to run the producer

    _consumers = []
    for cid in range(nconsumer):        # create each consumer and run it
        _consumer = _get_consumer(target=consumer, name=consumer_name,
                                  args=consumer_args, kwargs=consumer_kwargs,
                                  nproducer=nproducer,
                                  iqueue=producer_queue,
                                  oqueue=consumer_queue)
        _consumer.start()
        _consumers.append(_consumer)

    # create reporter and run it
    _reporter = _get_reporter(target=reporter, name=reporter_name,
                              args=reporter_args, kwargs=reporter_kwargs,
                              nconsumer=nconsumer, iqueue=consumer_queue)
    _reporter.start()

    _producer.join()
    for _consumer in _consumers:
        _consumer.join()
    _reporter.join()


