#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Kyle A. Barlow, Shane O'Connor
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

"""\
The script kicks off the preminmization step of the benchmark run using the minimize_with_cst application from the
Rosetta suite. The command lines used herein are intended to reproduce the protocol from row 16 of the original paper by Kellogg et al.:

 Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein
 structure and stability. 2011. Proteins. 79(3):830-8. doi: 10.1002/prot.22921.

Usage:
    run_preminimization.py [options]...

Options:

    -d --dataset DATASET
        A filepath to the input dataset in JSON format. [default: ../../input/json/kellogg.json]

    -o --output_directory OUTPUT_DIR
        The path where output data will be created. Output will be created inside a time-stamped subfolder of this directory. [default: ./job_output]

    --run_identifier RUN_ID
        A suffix used to name the output directory.

    --test
        When this option is set, a shorter version of the benchmark will run with fewer input structures, less fewer DDG experiments, and fewer generated structures. This should be used to test the scripts but not for analysis.

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
import rosetta.parse_settings
from rosetta.write_run_file import process as write_run_file
from analysis.libraries import docopt
from analysis.stats import read_file, write_file
try:
    import json
except:
    import simplejson as json

from rosetta.pdb import PDB, create_mutfile
from rosetta.basics import Mutation


task_subfolder = 'preminimization'
mutfiles_subfolder = 'mutfiles'
generated_scriptname = 'preminimization_step'

def create_input_files(pdb_dir_path, pdb_data_dir, mutfile_data_dir, keypair, dataset_cases, skip_if_exists = False):
    '''Create the stripped PDB files and the mutfiles for the DDG step. Mutfiles are created at this point as we need the
    original PDB to generate the residue mapping.
    '''

    # Read PDB
    pdb_id = keypair[0]
    chain = keypair[1]
    pdb = PDB.from_filepath(pdb_dir_path)
    stripped_pdb_path = os.path.join(pdb_data_dir, '%s_%s.pdb' % (pdb_id, chain))

    # Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
    chains = [chain]
    pdb.stripForDDG(chains, numberOfModels = 1)

    # Check to make sure that we haven't stripped all the ATOM lines
    if not [line for line in pdb.lines if line[0:4] == "ATOM"]:
        raise Exception("No ATOM lines remain in the stripped PDB file %s." % stripped_pdb_path)

    # Check to make sure that CSE and MSE are not present in the PDB
    unhandled_residues = pdb.CheckForPresenceOf(["CSE", "MSE"])
    if unhandled_residues:
        raise Exception("Found residues [%s] in the stripped PDB file %s. These should be changed to run this job under Rosetta." % (', '.join(unhandled_residues), stripped_pdb_path))

    # Turn the lines array back into a valid PDB file
    if not(skip_if_exists) or not(os.path.exists(stripped_pdb_path)):
        write_file(stripped_pdb_path, '\n'.join(pdb.lines))

    # Check that the mutated positions exist and that the wild-type matches the PDB
    for dataset_case in dataset_cases:
        # todo: Hack. This should be removed when PDB homologs are dealt with properly.
        mutation_objects = []
        assert(dataset_case['PDBFileID'] == pdb_id)
        for mutation in dataset_case['Mutations']:
            #if experimentPDB_ID == "1AJ3" and predictionPDB_ID == "1U5P":
            #    assert(int(mutation['ResidueID']) < 1000)
            #    mutation['ResidueID'] = str(int(mutation['ResidueID']) + 1762)
            mutation_objects.append(Mutation(mutation['WildTypeAA'], mutation['ResidueID'], mutation['MutantAA'], mutation['Chain']))

        pdb.validate_mutations(mutation_objects)

        # Post stripping checks
        # Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
        # then check again that the mutated positions exist and that the wild-type matches the PDB

        remappedMutations = pdb.remapMutations(mutation_objects, pdb_id)
        mutfile = create_mutfile(pdb, remappedMutations)
        mutfilename = os.path.join(mutfile_data_dir, '%d.mutfile' % (dataset_case['RecordID']))
        if os.path.exists(mutfilename):
            raise Exception('%s already exists. Check that the RecordIDs in the JSON file are all unique.' % mutfilename)
        write_file(os.path.join(mutfile_data_dir, '%d.mutfile' % (dataset_case['RecordID'])), mutfile)
    return stripped_pdb_path


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
    settings = rosetta.parse_settings.get_dict()

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

    # Set the job directory name
    job_name = '%s_%s_ddg_monomer_16' % (time.strftime("%y-%m-%d-%H-%M"), getpass.getuser())
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
    output_dir = os.path.join(root_output_directory, job_name) # The root directory for the protocol run
    task_dir = os.path.join(output_dir, task_subfolder) # The root directory for preminization section of the protocol
    output_data_dir = os.path.join(output_dir, 'data')
    pdb_data_dir = os.path.join(output_data_dir, 'input_pdbs')
    mutfile_data_dir = os.path.join(output_data_dir, mutfiles_subfolder)
    for jobdir in [output_dir, task_dir, output_data_dir, pdb_data_dir, mutfile_data_dir]:
        try: os.mkdir(jobdir)
        except: pass

    # Make a copy the dataset so that it can be automatically used by the following steps
    shutil.copy(dataset_filepath, os.path.join(output_dir, 'dataset.json'))

    # Count the number of datapoints per PDB chain
    count_by_pdb_chain = {}
    dataset_cases_by_pdb_chain = {}
    job_dict = {}
    for ddg_case in dataset_cases:
        chains = set([r['Chain'] for r in ddg_case['Mutations']])
        assert(len(chains) == 1)
        chain = chains.pop()
        pdb_id = ddg_case['PDBFileID']
        keypair = (pdb_id, chain)
        count_by_pdb_chain[keypair] = count_by_pdb_chain.get(keypair, 0)
        count_by_pdb_chain[keypair] += 1
        dataset_cases_by_pdb_chain[keypair] = dataset_cases_by_pdb_chain.get(keypair, [])
        dataset_cases_by_pdb_chain[keypair].append(ddg_case)

    # Create the list of PDB IDs and chains for the dataset
    print('')
    if arguments['--test']:
        pdb_monomers = []
        print('Creating test run input...')
        num_cases = 0
        for keypair, v in sorted(count_by_pdb_chain.iteritems(), key=lambda x:-x[1]):
            if v <= 10:
                pdb_monomers.append(keypair)
                num_cases += v
                if num_cases >= 20:
                    break
    else:
        pdb_monomers = sorted(count_by_pdb_chain.keys())

    # Ensure all the input PDB files exist
    for keypair in pdb_monomers:
        pdb_path = os.path.join(input_pdb_dir_path, '%s.pdb' % keypair[0])
        if not os.path.exists(pdb_path):
            raise Exception('Error: The file %s is missing.' % pdb_path)

    # Write job dict and setup self-contained data directory
    print('Copying PDB files:')
    job_dict = {}

    for keypair in pdb_monomers:
        sys.stdout.write('.')
        sys.stdout.flush()
        sub_dict = {}

        # Create stripped PDB file
        pdb_dir_path = os.path.join(input_pdb_dir_path, '%s.pdb' % keypair[0])
        #new_pdb_path = os.path.join(pdb_data_dir, os.path.basename(pdb_path))
        #if not os.path.isfile(new_pdb_path):
        #    # Note: We are assuming that the existing file has not been modified to save time copying the files
        #    shutil.copy(pdb_path, new_pdb_path)
        stripped_pdb_path = create_input_files(pdb_dir_path, pdb_data_dir, mutfile_data_dir, keypair, dataset_cases_by_pdb_chain[keypair])
        pdb_relpath = os.path.relpath(stripped_pdb_path, output_dir)

        # Set up --in:file:l parameter
        sub_dict['input_file_list'] = [pdb_relpath]

        job_dict[os.path.join(task_subfolder, '_'.join(keypair))] = sub_dict
    sys.stdout.write('\n')

    with open(os.path.join(output_data_dir, 'job_dict.pickle'), 'w') as f:
        pickle.dump(job_dict, f)

    settings['numjobs'] = '%d' % len(pdb_monomers)
    settings['mem_free'] = '3.0G'
    settings['scriptname'] = generated_scriptname
    settings['appname'] = 'minimize_with_cst'
    settings['rosetta_args_list'] =  [
        '-in:file:fullatom', '-ignore_unrecognized_res',
        '-fa_max_dis', '9.0', '-ddg::harmonic_ca_tether', '0.5',
        '-ddg::constraint_weight', '1.0',
        '-ddg::out_pdb_prefix', 'min_cst_0.5',
        '-ddg::sc_min_only', 'false'
    ]
    settings['output_dir'] = output_dir

    write_run_file(settings)
    job_path = os.path.abspath(output_dir)
    print('''Job files written to directory: %s. To launch this job:
    cd %s
    python %s''' % (job_path, job_path, generated_scriptname))


