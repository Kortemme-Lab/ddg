#!/usr/bin/env python2
# This work is licensed under the terms of the MIT license. See LICENSE for the full text.

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
