#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Shane O'Connor
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# This should be moved into the tools repository.


class ReportingObject(object):
    '''A simple class to allow stdout suppression.'''

    def __init__(self, silent = False):
        self.silent = silent


    def log(self, str, fn = None):
        if not self.silent:
            if fn:
                fn(str)
            else:
                print(str)


    def glog(self, *args, **kwargs):
        # A more generic logging function accepting args and kwargs
        if not self.silent:
            if 'fn' in kwargs:
                fn = kwargs['fn']
                del kwargs['fn']
                fn(*args, **kwargs)
                kwargs['fn'] = fn
            else:
                if args and kwargs: print(args, kwargs)
                elif kwargs: print(kwargs)
                else: print(args)