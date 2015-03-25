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

import os
import sys
import inspect
import json

path_to_this_module = os.path.abspath( os.path.dirname( inspect.getsourcefile(sys.modules[__name__]) ) )
settings_json = 'settings.json'

def get_dict():
    settings_file = os.path.join(path_to_this_module, settings_json)

    # Read the settings file
    if not os.path.exists(settings_file):
        raise Exception('The settings file settings.json does not exist (expected location: "%s"). Please create this file appropriately; see %s/settings.json.example contains an example configuration file which can be used as a template.' % (os.path.abspath(settings_file), os.path.split(os.path.abspath(settings_file))[0]))
    try:
        with open(settings_file, 'r') as f:
            settings = json.load(f)
    except Exception, e:
        raise Exception('An error occurred parsing the settings file: %s.\n%s' % (settings_json, str(e)))

    required_settings = [
        "local_rosetta_installation_path",
        "cluster_rosetta_installation_path",
        "rosetta_binary_type",
    ]

    for setting in required_settings:
        if setting not in settings:
            raise Exception('The setting argument "%s" is required in %s; please add' % (setting, settings_json) ) 

    settings['cluster_rosetta_bin'] = os.path.join(settings['cluster_rosetta_installation_path'], 'source', 'bin')
    settings['cluster_rosetta_db_dir'] = os.path.join(settings['cluster_rosetta_installation_path'], 'database')

    settings['local_rosetta_bin'] = os.path.join(settings['local_rosetta_installation_path'], 'source', 'bin')
    settings['local_rosetta_db_dir'] = os.path.join(settings['local_rosetta_installation_path'], 'database')

    return settings
