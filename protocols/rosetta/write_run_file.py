#!/usr/bin/python

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
            
    # Handle LD paths
    if 'extra_ld_path' in data_dict:
        data_dict['add_extra_ld_path'] = 'True'
    else:
        data_dict['add_extra_ld_path'] = 'False'
        data_dict['extra_ld_path'] = ''

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

