#!/usr/bin/python

import os
import sys

template_file = 'template.py'

def process(data_dict):
    required_arguments = [
        'numjobs', 'scriptname', 'mem_free',
        'cluster_rosetta_bin', 'cluster_rosetta_db',
        'local_rosetta_bin', 'local_rosetta_db',
        'appname', 'rosetta_args_list',
        'output_dir',
    ]

    # Arguments that cannot be in data_dict
    unrequired_arguments = [
        'add_extra_ld_path',
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

    if not os.path.isdir(data_dict['output_dir']):
        os.makedirs(data_dict['output_dir'])

    formatted_data_dict = {}
    for arg in data_dict:
        new_arg = '#$#%s#$#' % arg
        value = data_dict[arg]

        formatted_data_dict[new_arg] = data_dict[arg]

    new_lines = []
    with open(template_file, 'r') as f:
        for line in f:
            for arg in formatted_data_dict:
                line = line.replace(arg, formatted_data_dict[arg])
            new_lines.append(line)

    with open(os.path.join(data_dict['output_dir'], '%s.py' % data_dict['scriptname']), 'w') as f:
        for line in new_lines:
            f.write(line)

