#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Kyle A. Barlow
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

import os
import sys
import inspect
import math

path_to_this_module = os.path.abspath( os.path.dirname( inspect.getsourcefile(sys.modules[__name__]) ) )

template_file = os.path.join(path_to_this_module, 'template.py')

def process(data_dict):
    required_arguments = [
        'numjobs', 'scriptname', 'mem_free',
        'cluster_rosetta_bin',
        'local_rosetta_bin',
        'appname', 'rosetta_args_list',
        'output_dir',
        'rosetta_binary_type',
    ]

    # Arguments that cannot be in data_dict
    unrequired_arguments = [
        'add_extra_ld_path',
        'numclusterjobs',
    ]

    for arg in required_arguments:
        if arg not in data_dict:
            print 'ERROR: Data dictionary missing argument', arg
            sys.exit(1)

    for arg in unrequired_arguments:
        if arg in data_dict:
            print 'ERROR: Data dictionary cannot contain argument', arg
            sys.exit(1)

    # Complete appname
    data_dict['fullappname'] = data_dict['appname'] + data_dict['rosetta_binary_type']
            
    # Handle LD paths
    if 'extra_ld_path' in data_dict:
        data_dict['add_extra_ld_path'] = 'True'
    else:
        data_dict['add_extra_ld_path'] = 'False'
        data_dict['extra_ld_path'] = ''

    # Handle if general rosetta args are a list instead of a string
    if not isinstance(data_dict['rosetta_args_list'], basestring):
        rosetta_args = ""
        if len(data_dict['rosetta_args_list']) == 1:
            rosetta_args += "'%s'" % data_dict['rosetta_args_list']
        elif len(data_dict['rosetta_args_list']) > 1:
            for arg in data_dict['rosetta_args_list'][:-1]:
                rosetta_args += "'%s', " % arg
            rosetta_args += "'%s'" % data_dict['rosetta_args_list'][-1]
        data_dict['rosetta_args_list'] = rosetta_args

    # Handle other options
    if 'cluster_rosetta_db' not in data_dict:
        data_dict['cluster_rosetta_db'] = ''

    if 'local_rosetta_db' not in data_dict:
        data_dict['local_rosetta_db'] = ''

    if 'tasks_per_process' not in data_dict:
        data_dict['tasks_per_process'] = 1
        data_dict['numclusterjobs'] = data_dict['numjobs']
    elif data_dict['tasks_per_process'] > 1:
        data_dict['numclusterjobs'] = int(
            math.ceil( float(data_dict['numjobs']) / float(data_dict['tasks_per_process']) )
        )

    if not os.path.isdir(data_dict['output_dir']):
        os.makedirs(data_dict['output_dir'])

    formatted_data_dict = {}
    for arg in data_dict:
        new_arg = '#$#%s#$#' % arg
        value = data_dict[arg]

        formatted_data_dict[new_arg] = str(data_dict[arg])

    new_lines = []
    with open(template_file, 'r') as f:
        for line in f:
            for arg in formatted_data_dict:
                line = line.replace(arg, formatted_data_dict[arg])
            new_lines.append(line)

    with open(os.path.join(data_dict['output_dir'], '%s.py' % data_dict['scriptname']), 'w') as f:
        for line in new_lines:
            f.write(line)

