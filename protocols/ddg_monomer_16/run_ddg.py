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
The script kicks off the ddg step of the benchmark run using the ddg_monomer application from the
Rosetta suite. The command lines used herein are intended to reproduce the protocol from row 16 of the original paper by Kellogg et al.:

 Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein
 structure and stability. 2011. Proteins. 79(3):830-8. doi: 10.1002/prot.22921.

Usage:
    run_ddg.py [options]...

Options:

    -o --output_directory OUTPUT_DIR
        The path to a directory previously generated from the run_preminimization script. This defaults to the most recent directory in job_output, if this exists.

    -n --num_struct NUM_STRUCT
        This specifies the number of wildtype/mutant structures generated. If this is used with --test then the --test value for this option takes priority. [default: 50]

    --force
        When this option is set, the most recent directory in job_output, if it exists, will be used without prompting the user.

    --test
        When this option is set, a shorter version of the benchmark will run with fewer input structures, less fewer DDG experiments, and fewer generated structures. This should be used to test the scripts but not for analysis.

    --talaris2014
        When this option is set, the talaris2014 score function will be used rather than the default score function. Warning: This option may break when talaris2014 becomes the default Rosetta score function.

Authors:
    Kyle Barlow
    Shane O'Connor
"""

import sys
import os
import shutil
import glob
import cPickle as pickle
try: import json
except: import simplejson as json

from libraries import docopt

import rosetta.parse_settings
from rosetta.write_run_file import process as write_run_file

from klab.fsio.fs import read_file, write_file
from klab.tui.utils import prompt_yn
from run_preminimization import task_subfolder as preminimization_task_subfolder, mutfiles_subfolder


task_subfolder = 'ddg'
generated_scriptname = 'ddg_step'


def create_constraints_file(preminimization_log, outfile_path):
    '''This does the work of the convert_to_cst_file.sh script in the Rosetta repository.'''
    constraints = []
    contents = read_file(preminimization_log)
    for line in contents.split('\n'):
        if line.startswith("c-alpha"):
            line = line.split()
            constraints.append("AtomPair CA %s CA %s HARMONIC %s %s" % (line[5], line[7], line[9], line[12]))
    write_file(outfile_path, '\n'.join(constraints))
    return outfile_path


def create_constraints_files(preminimized_pdb_data_dir, constraints_data_dir):
    constraints_files = {}
    preminimized_structures = {}
    for pdir in os.listdir(preminimized_pdb_data_dir):
        if len(pdir) == 6 and pdir[4] == '_':
            pdb_id = pdir.split('_')[0]
            chain_id = pdir.split('_')[1]
            pcase_path = os.path.join(preminimized_pdb_data_dir, pdir)
            output_files = os.listdir(pcase_path)
            output_structure = 'min_cst_0.5.%s_0001.pdb.gz' % pdir
            if 'rosetta.out.gz' not in output_files:
                raise Exception('The expected output file rosetta.out.gz was not found in %s.' % pcase_path)
            if (output_structure) not in output_files:
                raise Exception('The expected preminimized structure %s was not found in %s.' % (output_structure, pcase_path))
            constraints_files[(pdb_id, chain_id)] = create_constraints_file(os.path.join(pcase_path, 'rosetta.out.gz'), os.path.join(constraints_data_dir, '%s_%s.cst' % (pdb_id, chain_id)))
            preminimized_structures[(pdb_id, chain_id)] = os.path.join(pcase_path, output_structure)
    return constraints_files, preminimized_structures


if __name__ == '__main__':
    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)

    # Determine the output directory
    output_dir = None
    if arguments.get('--output_directory'):
        output_dir = arguments['--output_directory'][0]
        if not(os.path.exists(output_dir)):
            raise Exception('The directory %s does not exist.' % output_dir)
    else:
        output_dir = os.path.abspath('job_output')
        if os.path.exists(output_dir):
            existing_dirs = [os.path.join(output_dir, d) for d in os.listdir(output_dir) if d.find('ddg_monomer_16')!=-1 and os.path.isdir(os.path.join(output_dir, d))]
            most_recent_directory = sorted(existing_dirs)[-1]
            if most_recent_directory:
                answer = None
                if arguments.get('--force'):
                    answer = True
                    print('\nRunning the ddg_monomer step in %s.' % most_recent_directory)
                else:
                    answer = prompt_yn('\nNo output path was specified. Use %s (y/n)?' % most_recent_directory)
                if not answer:
                    print('No output path was specified. Exiting.\n')
                    sys.exit(1)
                output_dir = most_recent_directory
            else:
                print('No preminimization output could be found in the job_output directory. Exiting.\n')
                sys.exit(1)

    # Read the settings file
    settings = rosetta.parse_settings.get_dict()

    # Set the job output directories
    output_data_dir = os.path.join(output_dir, 'data')
    mutfile_data_dir = os.path.join(output_data_dir, mutfiles_subfolder)
    preminimized_pdb_data_dir = os.path.join(output_dir, preminimization_task_subfolder)
    constraints_data_dir = os.path.join(output_data_dir, 'constraints')
    try: os.mkdir(constraints_data_dir)
    except: pass
    for p in [output_data_dir, mutfile_data_dir, preminimized_pdb_data_dir, constraints_data_dir]:
        if not os.path.exists(p):
            raise Exception('The folder %s was expected to exist after the preminimization step but could not be found.' % p)

    # Read in the dataset file
    try:
        dataset_filepath = os.path.join(output_dir, 'dataset.json')
        dataset = json.loads(read_file(dataset_filepath))
        dataset_cases = dataset['data']
    except Exception, e:
        raise Exception('An error occurred parsing the JSON file: %s..' % str(e))

    # Run all cases with associated mutfiles
    mutfiles = glob.glob(os.path.join(mutfile_data_dir, '*.mutfile'))
    existing_case_ids = []
    for m in mutfiles:
        try:
            filename = os.path.split(m)[1]
            assert(filename.endswith('.mutfile'))
            record_id = int(filename[:-8])
            existing_case_ids.append(record_id)
        except:
            raise Exception('The file %s was expected to have a name like [record_id].mutfile e.g. 1245.mutfile.' % m)

    count_by_pdb = {}
    job_dict = {}
    dataset_cases_by_id = {}
    for ddg_case in dataset_cases:
        dataset_cases_by_id[ddg_case['RecordID']] = ddg_case
    for existing_case_id in existing_case_ids:
        if existing_case_id not in dataset_cases_by_id:
            raise Exception('The dataset case corresponding to %d.mutfile could not be found in %s.' % (existing_case_id, dataset_filepath))

    # Write job dict and setup self-contained data directory
    extra_s = ''
    if arguments['--talaris2014']:
        extra_s = ' (using talaris2014)'
    print('Creating constraint files...%s' % extra_s)
    constraints_files, preminimized_structures = create_constraints_files(preminimized_pdb_data_dir, constraints_data_dir)

    number_of_structural_pairs = arguments['--num_struct'][0]
    if arguments['--test']:
        number_of_structural_pairs = 2 # only create two wildtype/mutant pairs in test mode

    for existing_case_id in existing_case_ids:
        dataset_case = dataset_cases_by_id[existing_case_id]
        pdb_id = dataset_case['PDBFileID']
        chains = set([m['Chain'] for m in dataset_case['Mutations']])
        if not len(chains) == 1:
            raise Exception('This script is only intended for monomeric structures but the set of mutations in case %d of %s uses more than one chain.' % (existing_case_id, dataset_filepath))
        chain_id = chains.pop()

        sys.stdout.write('.')
        sys.stdout.flush()

        sub_dict = {}
        constraints_file = constraints_files.get((pdb_id, chain_id))
        preminimized_structure = preminimized_structures.get((pdb_id, chain_id))
        mutfile = os.path.join(mutfile_data_dir, '%d.mutfile' % existing_case_id)
        if not constraints_file:
            raise Exception('Could not determine the constraints file for %s, chain %s.' % (pdb_id, chain_id))
        if not preminimized_structure:
            raise Exception('Could not determine the preminimized structure file for %s, chain %s.' % (pdb_id, chain_id))
        if not os.path.exists(mutfile):
            raise Exception('Could not locate the mutfile %s for dataset record %d.' % (mutfile, existing_case_id))
        sub_dict['-constraints::cst_file'] = os.path.relpath(constraints_file, output_dir)
        sub_dict['-in:file:s'] = os.path.relpath(preminimized_structure, output_dir)
        sub_dict['-ddg::mut_file'] = os.path.relpath(mutfile, output_dir)

        job_dict[os.path.join(task_subfolder, str(existing_case_id))] = sub_dict
    sys.stdout.write('\n')

    # Keep a copy of the preminimization step pickle for debugging
    pickle_file = os.path.join(output_data_dir, 'job_dict.pickle')
    if os.path.exists(pickle_file):
        existing_job_keys = pickle.load(open(pickle_file, 'r')).keys()
        for k in existing_job_keys:
            if k.startswith(preminimization_task_subfolder):
                shutil.copy(pickle_file, os.path.join(output_data_dir, '%s_step_dict.pickle' % preminimization_task_subfolder))
                break

    with open(pickle_file, 'w') as f:
        pickle.dump(job_dict, f)

    settings['numjobs'] = '%d' % len(existing_case_ids)
    settings['mem_free'] = '5.0G'
    settings['scriptname'] = generated_scriptname
    settings['appname'] = 'ddg_monomer'
    settings['rosetta_args_list'] = [
        '-in:file:fullatom', '-ignore_unrecognized_res', '-fa_max_dis', '9.0',
        '-ddg::dump_pdbs' ,'true', '-ddg::suppress_checkpointing' ,'true',
        '-ddg:weight_file' ,'soft_rep_design' ,'-ddg::iterations' ,str(number_of_structural_pairs),
        '-ddg::local_opt_only' ,'false' ,'-ddg::min_cst' ,'true',
        '-ddg::mean' ,'false' ,'-ddg::min', 'true',
        '-ddg::sc_min_only' ,'false',
        '-ddg::ramp_repulsive', 'true'
    ]
    if arguments['--talaris2014']:
        settings['rosetta_args_list'].extend(['-talaris2014', 'true'])
    settings['output_dir'] = output_dir

    write_run_file(settings)
    job_path = os.path.abspath(output_dir)
    print('''Job files written to directory: %s.\n\nTo launch this job locally (this will take some time):
    cd %s
    python %s.py

It is recommended to run this on an SGE cluster in which case use these commands instead:
    cd %s
    qsub %s.py\n''' % (job_path, job_path, generated_scriptname, job_path, generated_scriptname))


