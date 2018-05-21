#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : spider
# AUTHOR  : codeunsolved@gmail.com
# CREATED : February 15 2017
# VERSION : v0.0.1
# UPDAYE  : [v0.0.1] May 18 2018
# 1. optimize :AsyncSpider:'s INPUT and structure;

import time
from collections import defaultdict

from tornado import httpclient, gen, ioloop, queues

from .base import color_term


class AsyncSpider(object):
    """AsyncSpider.

    :param unfetched: dict{label: url}
    :param handle_response: func(response, label, url)

    :Refer:
    [example - a concurrent web spider](http://www.tornadoweb.org/en/stable/guide/queues.html?highlight=spider)
    """

    def __init__(self, unfetched,
                 handle_response,
                 concurrency=3, max_trial=10):
        self.unfetched = unfetched
        self.labels = sorted(self.unfetched.keys())
        self.handle_response = handle_response

        self.q = queues.Queue()
        self.fetching, self.fetched = set(), set()
        
        self.concurrency = concurrency
        self.max_trial = max_trial

        self.failed_trials = defaultdict(int)
        self.abandoned_trials = []

        self.start = None

        self.check_handle_response()
        [self.q.put(x) for x in self.labels]

    def check_handle_response(self):
        if not hasattr(self.handle_response, '__call__'):
            color_term("`handle_response()` is invalid!", 'ERR', exit(1))

    @gen.coroutine
    def get(self, label, url):
        try:
            resp = yield httpclient.AsyncHTTPClient().fetch(url)
        except Exception as e:
            color_term("Exception: {} for '{}'({})".format(e, label, url), 'ERR')
            raise gen.Return([])
        else:
            raise gen.Return(resp)

    @gen.coroutine
    def query(self):
        @gen.coroutine
        def fetch(i):
            label = yield self.q.get()
            url = self.unfetched[label]
            try:
                color_term("- [{}/{}]@{} {}".format(self.labels.index(label) + 1, len(self.labels), i, label))
                if label in self.fetching:
                    color_term("{} is fetching!".format(label), 'ERR')
                else:
                    self.fetching.add(label)
                    response = yield self.get(label, url)
                    if response:
                        self.handle_response(response, label, url)
                        self.fetched.add(label)
                    else:
                        self.failed_trials[label] += 1
                        if self.failed_trials[label] < self.max_trial:
                            color_term("FAILED! put {} back to the queue".format(label), 'WRN')
                            self.fetching.remove(label)
                            self.q.put(label)
                        else:
                            color_term("{} have reached max trails({})".format(label, self.max_trial), 'WRN')
                            self.fetched.add(label)
                            self.abandoned_trials.append(label)
            finally:
                self.q.task_done()

        @gen.coroutine
        def worker(i):
            while True:
                yield fetch(i)

        for _ in range(self.concurrency):
            worker(_)

        yield self.q.join()
        assert self.fetching == self.fetched
        print('------------------')

        if self.unfetched:
            color_term("Done in {:.0f} min {:.2f} sec, fetched {} URLs".format(
                *(divmod(time.time() - self.start, 60) + (len(self.unfetched),))), 'green')
        else:
            color_term("Nothing need to fetch")

        if self.abandoned_trials:
            print('------------------')
            color_term("Abandoned trials:", 'ERR')
            for x in self.abandoned_trials:
                print(x)

    def run(self):
        color_term("• start AsyncSpider")
        self.start = time.time()
        io_loop = ioloop.IOLoop.current()
        io_loop.run_sync(self.query)
