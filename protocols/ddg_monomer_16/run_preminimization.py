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

    --talaris2014
        When this option is set, the talaris2014 score function will be used rather than the default score function. Warning: This option may break when talaris2014 becomes the default Rosetta score function.

    --beta_july15
        When this option is set, the July 2015 beta score function will be used rather than the default score function. Warning: This option may break when this score function is removed.

    -p --parallel NUM_PROCESSORS
        If this argument is set then the job setup will use NUM_PROCESSORS which will speed this step up. Otherwise, a single processor will be used. This should run on both Unix and Windows machines.

    --maxp
        This is a special case of --parallel. If this argument is set then the job setup will use as many processors as are available on the machine.

Authors:
    Kyle Barlow
    Shane O'Connor
"""

import sys
import os
import shutil
import traceback
import time
import cPickle as pickle
import getpass

multiprocessing_module_available = True
try: import multiprocessing
except: multiprocessing_module_available = False
try: import json
except: import simplejson as json

from libraries import docopt

import rosetta.parse_settings
from rosetta.write_run_file import process as write_run_file

from klab.fs.fsio import read_file, write_file
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation
from klab.rosetta.input_files import Mutfile


task_subfolder = 'preminimization'
mutfiles_subfolder = 'mutfiles'
generated_scriptname = 'preminimization_step'
DEFAULT_NUMBER_OF_PROCESSORS_TO_USE = 1 # this can be overridden using the -p command


def create_input_files(job_dict, settings, pdb_dir_path, pdb_data_dir, mutfile_data_dir, keypair, dataset_cases, skip_if_exists = False):
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
    pdb.strip_to_chains(chains)
    pdb.strip_HETATMs()
    stripped_pdb = PDB('\n'.join(pdb.lines))

    # Check to make sure that we haven't stripped all the ATOM lines
    if not [line for line in stripped_pdb.lines if line[0:4] == "ATOM"]:
        raise Exception("No ATOM lines remain in the stripped PDB file %s." % stripped_pdb_path)

    # Assert that there are no empty sequences
    assert(sorted(stripped_pdb.atom_sequences.keys()) == sorted(chains))
    for chain_id, sequence in stripped_pdb.atom_sequences.iteritems():
        assert(len(sequence) > 0)

    # Check for CSE and MSE
    try:
        if 'CSE' in stripped_pdb.residue_types:
            raise Exception('This case contains a CSE residue which may (or may not) cause an issue with Rosetta depending on the version.')
        elif 'MSE' in stripped_pdb.residue_types:
            raise Exception('This case contains an MSE residue which may (or may not) cause an issue with Rosetta depending on the version.')
            # It looks like MSE (and CSE?) may now be handled - https://www.rosettacommons.org/content/pdb-files-rosetta-format
    except Exception, e:
        print('%s: %s, chain %s' % (str(e), str(stripped_pdb.pdb_id), chain))

    # Turn the lines array back into a valid PDB file
    if not(skip_if_exists) or not(os.path.exists(stripped_pdb_path)):
        write_file(stripped_pdb_path, '\n'.join(stripped_pdb.lines))

    # Create the mapping between PDB and Rosetta residue numbering
    # Note: In many Rosetta protocols, '-ignore_unrecognized_res' and '-ignore_zero_occupancy false' are used to allow
    # Rosetta to work with structures with missing data and non-canonicals. In those cases, we should supply both flags
    # in the string below. Since protocol 16 only uses '-ignore_unrecognized_res', we only use that flag below as otherwise
    # we could break the mapping.
    rosetta_scripts_bin = os.path.join(settings['local_rosetta_bin'], 'rosetta_scripts%s' % settings['rosetta_binary_type'])
    rosetta_database_path = settings['local_rosetta_db_dir']
    if not os.path.exists(rosetta_scripts_bin):
        raise Exception('The Rosetta scripts executable "{0}" could not be found. Please check your configuration file.'.format(rosetta_database_path))
    if not os.path.exists(rosetta_database_path):
        raise Exception('The path to the Rosetta database "{0}" could not be found. Please check your configuration file.'.format(rosetta_database_path))
    stripped_pdb.construct_pdb_to_rosetta_residue_map(rosetta_scripts_bin,rosetta_database_path, extra_command_flags = '-ignore_unrecognized_res')
    atom_to_rosetta_residue_map = stripped_pdb.get_atom_sequence_to_rosetta_json_map()
    rosetta_to_atom_residue_map = stripped_pdb.get_rosetta_sequence_to_atom_json_map()

    # Save the PDB <-> Rosetta residue mappings to disk
    write_file(os.path.join(pdb_data_dir, '%s_%s.rosetta2pdb.resmap.json' % (pdb_id, chain)), rosetta_to_atom_residue_map)
    write_file(os.path.join(pdb_data_dir, '%s_%s.pdb2rosetta.resmap.json' % (pdb_id, chain)), atom_to_rosetta_residue_map)

    # Assert that there are no empty sequences in the Rosetta-processed PDB file
    total_num_residues = 0
    d = json.loads(rosetta_to_atom_residue_map)
    for chain_id in chains:
        num_chain_residues = len([z for z in d.values() if z[0] == chain_id])
        total_num_residues += num_chain_residues
        assert(num_chain_residues > 0)

    # Check that the mutated positions exist and that the wild-type matches the PDB
    try:
        for dataset_case in dataset_cases:
            assert(dataset_case['PDBFileID'] == pdb_id)

            # Note: I removed a hack here for 1AJ3->1U5P mapping
            # The JSON file does not have the residue IDs in PDB format (5 characters including insertion code) so we need to repad them for the mapping to work
            pdb_mutations = [ChainMutation(mutation['WildTypeAA'], PDB.ResidueID2String(mutation['ResidueID']), mutation['MutantAA'], Chain = mutation['Chain']) for mutation in dataset_case['Mutations']]
            stripped_pdb.validate_mutations(pdb_mutations)

            # Map the PDB mutations to Rosetta numbering which is used by the mutfile format
            rosetta_mutations = stripped_pdb.map_pdb_residues_to_rosetta_residues(pdb_mutations)
            if (len(rosetta_mutations) != len(pdb_mutations)) or (None in set([m.ResidueID for m in rosetta_mutations])):
                raise Exception('An error occurred in the residue mapping code for DDG case: %s, %s' % (pdb_id, pdb_mutations))

            # Create the mutfile
            mutfile = Mutfile.from_mutagenesis(rosetta_mutations)
            mutfilename = os.path.join(mutfile_data_dir, '%d.mutfile' % (dataset_case['RecordID']))
            if os.path.exists(mutfilename):
                raise Exception('%s already exists. Check that the RecordIDs in the JSON file are all unique.' % mutfilename)
            write_file(os.path.join(mutfile_data_dir, '%d.mutfile' % (dataset_case['RecordID'])), str(mutfile))
    except Exception, e:
        print(str(e))
        print(traceback.format_exc())

    # Set up --in:file:l parameter
    pdb_relpath = os.path.relpath(stripped_pdb_path, settings['output_dir'])
    job_dict[os.path.join(task_subfolder, '_'.join(keypair))] = dict(input_file_list = [pdb_relpath])
    sys.stdout.write('.'); sys.stdout.flush()


def single_job_pack(args):
    print(args)
    return single_job(*args)


def use_multiple_processors(settings, pdb_monomers, input_pdb_dir_path, pdb_data_dir, mutfile_data_dir, dataset_cases_by_pdb_chain, num_processors):
    assert(multiprocessing_module_available)

    pool = multiprocessing.Pool(processes = num_processors)#[, initializer[, initargs]]])
    m = multiprocessing.Manager()
    job_dict = m.dict()

    pool_jobs = []
    for keypair in pdb_monomers:
        pdb_dir_path = os.path.join(input_pdb_dir_path, '%s.pdb' % keypair[0])
        pool_jobs.append(pool.apply_async(create_input_files, (job_dict, settings, pdb_dir_path, pdb_data_dir, mutfile_data_dir, keypair, dataset_cases_by_pdb_chain[keypair])))
    pool.close()
    pool.join()
    sys.stdout.write('\n')
    return job_dict._getvalue()


def use_single_processor(settings, pdb_monomers, input_pdb_dir_path, pdb_data_dir, mutfile_data_dir, dataset_cases_by_pdb_chain):
    job_dict = {}
    for keypair in pdb_monomers:
        pdb_dir_path = os.path.join(input_pdb_dir_path, '%s.pdb' % keypair[0])
        create_input_files(job_dict, settings, pdb_dir_path, pdb_data_dir, mutfile_data_dir, keypair, dataset_cases_by_pdb_chain[keypair])
    sys.stdout.write('\n')
    return job_dict


if __name__ == '__main__':
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

    # Read in any parallel processing options
    num_system_processors = 1
    if multiprocessing_module_available:
        num_system_processors = multiprocessing.cpu_count()

    if arguments.get('--maxp'):
        num_processors = num_system_processors
    else:
        num_processors = min(DEFAULT_NUMBER_OF_PROCESSORS_TO_USE, num_system_processors)
        if arguments.get('--parallel'):
            valid_options = [int(x) for x in arguments['--parallel'] if x.isdigit()]
            if not valid_options:
                raise Exception('None of the arguments to --parallel are valid. The argument must be an integer between 1 and the number of processors (%d).' % num_system_processors)
            else:
                num_processors = max(valid_options)
        else:
            # If the user has not specified the number of processors, only one is selected, and more exist then let them know that this process may run faster
            if num_processors == 1 and num_system_processors > 1:
                print('The setup is configured to use one processor but this machine has %d processors. The --parallel or --maxp options may make this setup run faster.' % num_system_processors)
    if 1 > num_processors or num_processors > num_system_processors:
        raise Exception('The number of processors must be an integer between 1 and %d.' % num_system_processors)

    # Read the dataset from disk
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
    settings['output_dir'] = output_dir
    try:
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
        extra_s = ''
        if arguments['--talaris2014']:
            extra_s = ' (using talaris2014)'
        if arguments['--beta_july15']:
            assert(not(extra_s))
            extra_s = ' (using beta_july15)'
        print('Creating benchmark input:%s' % extra_s)

        if num_processors == 1:
            job_dict = use_single_processor(settings, pdb_monomers, input_pdb_dir_path, pdb_data_dir, mutfile_data_dir, dataset_cases_by_pdb_chain)
        else:
            print('Setting up the preminimization data using %d processors.' % num_processors)
            job_dict = use_multiple_processors(settings, pdb_monomers, input_pdb_dir_path, pdb_data_dir, mutfile_data_dir, dataset_cases_by_pdb_chain, num_processors)

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
        if arguments['--talaris2014']:
            settings['rosetta_args_list'].extend(['-talaris2014', 'true'])
        elif arguments['--beta_july15']:
            settings['rosetta_args_list'].extend(['-beta_july15'])

        write_run_file(settings)
        job_path = os.path.abspath(output_dir)
        print('''Job files written to directory: %s.\n\nTo launch this job:
        cd %s
        python %s.py\n''' % (job_path, job_path, generated_scriptname))
    except Exception, e:
        print('\nAn exception occurred setting up the preminimization step: "%s".' % str(e))
        sys.stdout.write('Removing the directory %s: ' % output_dir)
        try:
            shutil.rmtree(output_dir)
            print('done.\n')
        except Exception, e2:
            print('failed.\n')
            print(str(e2))
            print(traceback.format_exc())



