#!/usr/bin/env python2
# This work is licensed under the terms of the MIT license. See LICENSE for the full text.

"""\
The script kicks off a benchmark run using the ddg_monomer application from the Rosetta suite. The command lines used
herein are intended to reproduce the protocol from row 16 of the original paper by Kellogg et al.:

 Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein
 structure and stability. 2011. Proteins. 79(3):830-8. doi: 10.1002/prot.22921.

Usage:
    run_benchmark.py [options]...

Options:

    -d DATASET --dataset DATASET
        A filepath to the input dataset in JSON format [default: ../../input/json/kellogg.json]

    -i RUN_ID --run_identifier RUN_ID
        A suffix used to name the output directory.

    -o OUTPUT_DIR --output_directory OUTPUT_DIR
        The path where output data will be created. Output will be created inside a time-stamped subfolder of this directory [default: ./job_output]

    -t --test
        When this option is set, a shorter version of the benchmark will run with limited sampling. This should be used to test the scripts but not for analysis.

Authors:
    Kyle Barlow
    Shane O'Connor
"""

import sys
import os
import re
import shutil
import time
import datetime
import inspect
import multiprocessing
import cPickle as pickle
import getpass
from rosetta.write_run_file import process as write_run_file
from analysis.libraries import docopt
from analysis.stats import read_file
try:
    import json
except:
    import simplejson as json


def make_resfile(resfile_path, mutation_datum, tanja_id):
    # Make resfile
    mutation_info_index = mutation_datum.tanja_id_list.index(tanja_id)
    chain = mutation_datum.chain_list[mutation_info_index]
    pdb_res = mutation_datum.pdb_res_list[mutation_info_index]
    insertion_code = mutation_datum.insertion_code_list[mutation_info_index]

    with open(resfile_path, 'w') as f:
        f.write('NATRO\nEX 1 EX 2 EX 3\nSTART\n')
        f.write('%d%s %s PIKAA A\n' % (pdb_res, insertion_code, chain))


if __name__ == '__main__':
    import pprint
    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)

    # Set the PDB input path
    input_pdb_dir_path = '../../input/pdbs'

    # Read the settings file
    if not os.path.exists('settings.json'):
        raise Exception('The settings file settings.json does not exist. Please create this file appropriately (see settings.json.example).')
    try:
        settings = json.loads(read_file('settings.json'))
    except Exception, e:
        raise Exception('An error occurred parsing the settings file settings.json: %s.' % str(e))
    local_rosetta_path = settings['local_rosetta_installation_path']
    cluster_rosetta_path = settings['cluster_rosetta_installation_path']
    rosetta_binary_type = settings['rosetta_binary_type']
    local_rosetta_bin_dir = os.path.join(local_rosetta_path, 'source', 'bin')
    local_rosetta_db_dir = os.path.join(local_rosetta_path, 'database')
    cluster_rosetta_bin_dir = os.path.join(cluster_rosetta_path, 'source', 'bin')
    cluster_rosetta_db_dir = os.path.join(cluster_rosetta_path, 'database')

    # Read in the dataset file
    dataset_filepath = arguments['--dataset'][0]
    dataset_filename = os.path.splitext(os.path.split(dataset_filepath)[1])[0]
    if not os.path.exists(dataset_filepath):
        raise Exception('The dataset file %s does not exist.' % dataset_filepath)

    try:
        dataset = json.loads(read_file(dataset_filepath))
        dataset_cases = dataset['data']
    except Exception, e:
        raise Exception('An error occurred parsing the JSON file: %s..' % str(e))

    # Read in the dataset file
    job_name = '%s_%s_ddg_monomer_16' % (time.strftime("%y-%m-%d-%H-%M"), getpass.getuser())
    job_name = '%s_%s_ddg_monomer_16' % (time.strftime("%y-%m-%d"), getpass.getuser())

    if arguments.get('--run_identifier'):
        job_name += '_' + arguments['--run_identifier'][0]

    # Set the root output directory
    root_output_directory = 'job_output'
    if arguments.get('--output_directory'):
        root_output_directory = arguments['--output_directory'][0]
    if not os.path.exists(root_output_directory):
        print('Creating directory %s:' % root_output_directory)
        os.makedirs(root_output_directory)

    # Set the job output directory
    output_dir = os.path.join(root_output_directory, job_name)
    output_data_dir = os.path.join(output_dir, 'data')
    pdb_data_dir = os.path.join(output_data_dir, 'input_pdbs')
    resfile_data_dir = os.path.join(output_data_dir, 'resfiles')
    try: os.mkdir(output_dir)
    except: pass
    try: os.mkdir(output_data_dir)
    except: pass
    try: os.mkdir(pdb_data_dir)
    except: pass
    try: os.mkdir(resfile_data_dir)
    except: pass





    '''
    {u'AggregateType': u'SingleValue',
  u'DDG': 0.956,
  u'ExperimentalDDGs': [{u'DDG': 0.956022944550669,
                         u'DDGType': u'DDG_H2O',
                         u'LocationOfValueInPublication': u'Table 1',
                         u'Publication': u'PMID:16042382',
                         u'Temperature': 25.0,
                         u'pH': 7.0}],
  u'Mutations': [{u'Chain': u'A',
                  u'DSSPExposure': 0.14375,
                  u'DSSPSimpleSSType': u'S',
                  u'DSSPType': u'E',
                  u'MutantAA': u'A',
                  u'ResidueID': u'95',
                  u'WildTypeAA': u'V'}],
  u'PDBFileID': u'5AZU',
  u'RecordID': 1206,
  u'_DataSetDDGID': 11915,
  u'_ExperimentID': 113616},
    '''

    count_by_pdb = {}
    job_dict = {}
    for ddg_case in dataset_cases:
        pdb_id = ddg_case['PDBFileID']
        count_by_pdb[pdb_id] = count_by_pdb.get(pdb_id, 0)
        count_by_pdb[pdb_id] += 1

    print('')
    if arguments['--test']:
        print('Creating test run input...')
        subset_of_pdb_paths = []
        num_cases = 0
        for k, v in sorted(count_by_pdb.iteritems(), key=lambda x:-x[1]):
            if v <= 10:
                subset_of_pdb_paths.append(os.path.join(input_pdb_dir_path, '%s.pdb' % k))
                num_cases += v
                if num_cases >= 20:
                    break
        pdb_paths = subset_of_pdb_paths
    else:
        pdb_paths = sorted(set([os.path.join(input_pdb_dir_path, '%s.pdb' % ddg_case['PDBFileID']) for ddg_case in dataset_cases]))

    for pdb_path in pdb_paths:
        if not os.path.exists(pdb_path):
            raise Exception('Error: The file %s is missing.' % pdb_path)

    # Write job dict and setup self-contained data directory
    print('Copying PDB files...')
    job_dict = {}

    for pdb_path in pdb_paths:
        sub_dict = {}
        new_pdb_path = os.path.join(pdb_data_dir, os.path.basename(pdb_path))
        if not os.path.isfile(new_pdb_path):
            # Note: We are assuming that the existing file has not been modified to save time copying the files
            shutil.copy(pdb_path, new_pdb_path)
        pdb_relpath = os.path.relpath(new_pdb_path, output_dir)
        sub_dict['input_file_list'] = [pdb_relpath]

        pdb_name = os.path.basename(pdb_path).split('.')[0]
        assert( len(pdb_name) == 4)
        job_dict[pdb_name.upper()] = sub_dict

    with open(os.path.join(output_data_dir, 'job_dict.pickle'), 'w') as f:
        pickle.dump(job_dict, f)

    args = {
        'numjobs' : '%d' % len(pdb_paths),
        'mem_free' : '3.0G',
        'scriptname' : 'cstmin_run',
        'cluster_rosetta_bin' : cluster_rosetta_bin_dir,
        'cluster_rosetta_db' : cluster_rosetta_db_dir,
        'local_rosetta_bin' : local_rosetta_bin_dir,
        'local_rosetta_db' : local_rosetta_db_dir,
        'appname' : 'minimize_with_cst%s' % rosetta_binary_type,
        'rosetta_args_list' :  "'-in:file:fullatom', '-ignore_unrecognized_res', '-fa_max_dis 9.0', '-ddg::harmonic_ca_tether 0.5', '-ddg::constraint_weight 1.0', '-ddg::out_pdb_prefix min_cst_0.5', '-ddg::sc_min_only false'",
        'output_dir' : output_dir,
    }

    write_run_file(args)
    print 'Job files written to directory: %s.\n' % os.path.abspath(output_dir)



'''
jobIDsvar=${jobIDs[$SGE_TASK_ID]}
constraintscst_filevar=${constraintscst_file[$SGE_TASK_ID]}
infilesvar=${infiles[$SGE_TASK_ID]}
mutfilevar=${mutfile[$SGE_TASK_ID]}

ddg_monomer.static.linuxgccrelease -in:file:s $infilesvar -ddg::mut_file $mutfilevar -constraints::cst_file $constraintscst_filevar
-database rosetta_database -ignore_unrecognized_res -in:file:fullatom -fa_max_dis 9.0 -ddg::dump_pdbs true
-ddg::suppress_checkpointing true -ddg:weight_file soft_rep_design -ddg::iterations 50 -ddg::local_opt_only false -ddg::min_cst true -ddg
::mean false -ddg::min true -ddg::sc_min_only false -ddg::ramp_repulsive true
'''
