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
The script compiles the analysis files for the benchmark run and then invokes the analysis script.

 Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein
 structure and stability. 2011. Proteins. 79(3):830-8. doi: 10.1002/prot.22921.

Usage:
    run_analysis.py <analysis_directory> [-b BENCHMARK_RUN_DIRECTORY...] [-n BENCHMARK_RUN_NAME...] [options]

Options:

    -b --benchmark_run_directory BENCHMARK_RUN_DIRECTORY
        The path to a directory previously generated from the run_preminimization script which contains a completed
        benchmark run. This defaults to the most recent directory in job_output, if this exists. This number of these
        arguments should equal the number of --benchmark_run_name arguments.

    -n --benchmark_run_name BENCHMARK_RUN_NAME
        This argument is used to name each benchmark run specified with -b with a human-readable name and is used as a prefix
        to name the generated output files. This number of these arguments should equal the number of --benchmark_run_directory arguments.

    --use_published_data
        When this option is set, the prediction data from the protocol in row 16 of the paper by Kellogg et al. is included
        as a benchmark run. This allows users to compare results with the original publication.

    --use_existing_benchmark_data
        When this option is set, the benchmark_data.json file is not regenerated. This saves time on subsequent calls to this analysis script but is disabled by default for safety.

    --force
        When this option is set, the most recent directory in job_output, if it exists, will be used without prompting the user.

    -X --do_not_generate_plots
        When this option is set, the analysis files are generated and the metrics are computed but graphical plots are not generated.

    --burial_cutoff CUTOFF
        The cutoff below which a residue is considered buried (value should be between 0-1.0). [default: 0.25]

    --scx_cutoff SCX
        The cutoff (absolute value) below which an experimental DDG is considered neutral. The cutoff is currently symmetric
        (neutrality is defined as the region -SCX kcal/mol to SCX kcal/mol). [default: 1.0]

    --scy_cutoff SCY
        The cutoff (absolute value) below which a predicted DDG is considered neutral. The cutoff is currently symmetric
        (neutrality is defined as the region -SCY energy units to SCY energy units). [default: 1.0]

    --take_lowest N
        When this option is set, the average of the N lowest-scoring (most stable) mutant and wildtypes structures are used to calculate the DDG value. [default: 3]

    --use_single_reported_value
        By default, the analysis takes the three lowest-scoring mutant structures and the three lowest-scoring wildtype structures to calculate the DDG value. This approach was taken by Kellogg et al. and should reduce stochastic noise. If this option is set, the single value reported by ddg_monomer is used instead. We do not recommend using this option.

    --include_derived_mutations
        Some datasets contain duplicated datapoints in the form of derived mutations e.g. the mutation 107L ->1L63 in the Kellogg set is the reverse of the 1L63 -> 107L mutation and measurement. Including derived mutations creates bias in the analysis so we remove them by default. To include these derived values in the analysis, use this flag.

Authors:
    Shane O'Connor
    Kyle Barlow
"""

import sys
import os
import shutil
import time
import datetime
import inspect
import multiprocessing
import re
import glob
import pprint
import getpass
import gzip
import numpy
import pprint
import subprocess
import traceback
import shlex
import copy
import pandas
from rosetta.write_run_file import process as write_run_file
from analysis.libraries import docopt
from analysis.libraries import colortext
from analysis.stats import read_file, read_file_lines, write_file, prompt_yn, fraction_correct, get_xy_dataset_statistics, plot, format_stats_for_printing, RInterface

from run_ddg import task_subfolder as ddg_task_subfolder
try:
    import json
except:
    import simplejson as json




###

# todo: This methods are now deprecated since we store this information in the dataframe. Update the callers and remove these


def determine_SL_class(record):
    '''Returns:
        - SL if all of the mutations are from smaller (by volume) residues to larger residues
        - LS if all of the mutations are from larger residues to smaller residues
        - XX if all of the mutations are to and from residues of the same size e.g. M -> L
        - O  if none of the cases above apply (this occurs when records with multiple mutations mix the types above).'''

    mutation_types = set()
    for m in record['Mutations']:
        wt_vol, mut_vol = amino_acid_details[m['WildTypeAA']]['van der Waals volume'], amino_acid_details[m['MutantAA']]['van der Waals volume']

        if wt_vol < mut_vol:
            mutation_types.add('SL')
        elif wt_vol > mut_vol:
            mutation_types.add('LS')
        else:
            mutation_types.add('XX')
    if len(mutation_types) == 1:
        return mutation_types.pop()
    return 'O'


def has_G_or_P(record):
    '''Returns True if any of the wildtype or mutant residues are glycine or proline and False otherwise.'''
    for m in record['Mutations']:
        if m['WildTypeAA'] == 'G' or m['MutantAA'] == 'G' or m['WildTypeAA'] == 'P' or m['MutantAA'] == 'P': # G is more likely
            return True
    return False
####
















# Regex used to extract data from Rosetta output
ddGMover_regex = re.compile("^protocols.moves.ddGMover:\s*mutate\s*.*?\s*wildtype_dG\s*is:\s*.*?and\s*mutant_dG\s*is:\s*.*?\s*ddG\s*is:\s*(.*)$")

# Colors used in the plots
plot_colors = dict(
    neon_green = '#39FF14',
    brown = '#774400',
    purple = '#440077',
    cornflower_blue = '#6495ED',
    firebrick = '#B22222',
)

# Human-readable descriptions for the volume breakdown
by_volume_descriptions = dict(
    SL = 'small-to-large mutations',
    LS = 'large-to-small mutations',
    XX = 'no change in volume',
)


# This module contains two classes.
# The BenchmarkManager class reads in the script options, extracts the score data from the benchmark runs, and controls
# the analysis.
# The BenchmarkRun class creates a dataframe containing the raw data used for analysis. This dataframe is then used to
# performs analysis for the particular run or between another runs.



class BenchmarkManager(object):
    '''This class is responsible for extract the data from benchmark runs and creating BenchmarkRun objects which can then be called to perform analysis.'''


    def __init__(self, arguments):
        '''Read the command-line options and extract the raw data from the benchmark runs.'''

        # Declare object variables
        self.benchmark_run_data = {}
        self.analysis_directory = None
        self.take_lowest = None
        self.burial_cutoff = None
        self.stability_classication_x_cutoff, self.stability_classication_y_cutoff = None, None

        # Parse command-line arguments and extract/load the benchmark input data
        self.parse_arguments(arguments)
        self.setup()


    def parse_arguments(self, arguments):
        '''Read command-line options.'''

        # Output directory. The path where the analysis files should be created.
        self.analysis_directory = os.path.abspath(arguments['<analysis_directory>'])

        # Plot-generation option
        self.generate_plots = not(arguments['--do_not_generate_plots'])

        # take-lowest option
        try:
            assert(arguments['--take_lowest'].isdigit())
            self.take_lowest = int(arguments['--take_lowest'])
            assert(self.take_lowest)
        except:
            raise colortext.Exception('take_lowest must be a positive non-zero integer: "%s" was passed.' % arguments['--take_lowest'])

        # burial_cutoff option. Residues with a DSSP exposure greater than this value are considered exposed
        try:
            self.burial_cutoff = float(arguments['--burial_cutoff'])
            if 0.0 > self.burial_cutoff or self.burial_cutoff > 1.0: raise Exception()
        except ValueError, e:
            raise colortext.Exception('burial_cutoff must be a float value.')
        except Exception, e:
            raise colortext.Exception('burial_cutoff must be between 0.0 and 1.0 (inclusive).')

        # stability classication cutoff values
        try:
            self.stability_classication_x_cutoff = abs(float(arguments['--scx_cutoff']))
            self.stability_classication_y_cutoff = abs(float(arguments['--scy_cutoff']))
        except ValueError, e:
            raise colortext.Exception('The stability classification cutoffs (--scx_cutoff, --scy_cutoff) must be float values.')


    def setup(self):
        '''Load in the data for all of the specified benchmark_run_directories. This code works on a number of assumptions.
           Preconditions which should be met by the benchmark setup scripts:
             - the dataset used to run the benchmark is contained in a file called dataset.json in the output directory (this will be created when the benchmark run is set up)
             - the ddg results are stored in a subfolder <ddg_task_subfolder> (defined by the run script)
           Postconditions set up by this script:
             - the extracted benchmark data (analysis_data) is:
                - stored in a file called benchmark_data.json in the benchmark run root directory; or
                - created in that location by this script.
           The loaded data is stored inside BenchmarkRun objects for analysis and cross-analysis. For convenience, the
           analyze method of this class performs all possible single and pair-wise analysis.
        '''

        # Determine the benchmark run directories and user-specified names
        benchmark_run_directories = arguments['--benchmark_run_directory']
        use_published_data = arguments['--use_published_data']
        benchmark_run_names = arguments['--benchmark_run_name']

        # Check for uniqueness
        if benchmark_run_directories and (len(benchmark_run_directories) != len(set(benchmark_run_directories))):
            raise colortext.Exception('The --benchmark_run_directory options must be unique.')
        if benchmark_run_names and (len(benchmark_run_names) != len(set(benchmark_run_names))):
            raise colortext.Exception('The --benchmark_run_name options must be unique.')

        # If benchmarks were chosen, make sure the directories exist and have corresponding names. Otherwise, use the most recent run.
        if benchmark_run_directories:
            if len(benchmark_run_directories) != len(benchmark_run_names):
                raise colortext.Exception('Each benchmark_run_directory argument (there are {0}) must have a corresponding benchmark_run_name argument (there are {1}).'.format(len(benchmark_run_directories), len(benchmark_run_names)))
            for benchmark_run_directory in benchmark_run_directories:
                if not(os.path.exists(benchmark_run_directory)):
                    raise colortext.Exception('The directory %s does not exist.' % benchmark_run_directory)
        else:
            benchmark_run_directory = os.path.abspath('job_output')
            if os.path.exists(benchmark_run_directory):
                most_recent_directory = None
                existing_dirs = sorted([os.path.join(benchmark_run_directory, d) for d in os.listdir(benchmark_run_directory) if d.find('ddg_monomer_16') != -1 and os.path.isdir(os.path.join(benchmark_run_directory, d))])
                if existing_dirs:
                    most_recent_directory = existing_dirs[-1]
                if most_recent_directory:
                    answer = None
                    if arguments.get('--force'):
                        answer = True
                        print('\nRunning analysis from the run in %s.' % most_recent_directory)
                    else:
                        answer = prompt_yn('\nNo benchmark run path was specified. Use %s (y/n)?' % most_recent_directory)
                    if not answer:
                        raise colortext.Exception('No benchmark run was specified. Exiting.\n')
                    benchmark_run_directories = [most_recent_directory]
                    if len(benchmark_run_names) == 0:
                        benchmark_run_names.append('')
                        colortext.warning('No benchmark_run_name argument was specified. Defaulting to a blank name.')
                    elif len(benchmark_run_names) > 1:
                        raise colortext.Exception('No benchmark_run_directory was specified (one recent directory was chosen automatically) but multiple --benchmark_run_name arguments were specified.')
        assert(len(benchmark_run_directories) == len(benchmark_run_names))

        # At least one benchmark run must be chosen (including the published benchmark run)
        if not use_published_data and len(benchmark_run_directories) == 0:
            raise colortext.Exception('No benchmark runs were specified and no output could be found in the job_output directory. Exiting.\n')

        # Load in the data for each benchmark
        for x in range(len(benchmark_run_directories)):

            benchmark_run_directory = benchmark_run_directories[x]
            benchmark_run_name = benchmark_run_names[x]
            colortext.message('Setting up the analysis data for {0} in {1}.'.format(benchmark_run_name, benchmark_run_directory))

            # Check that the job output directories exist
            output_data_dir = os.path.join(benchmark_run_directory, 'data')
            ddg_data_dir = os.path.join(benchmark_run_directory, ddg_task_subfolder)
            for p in [output_data_dir, ddg_data_dir]:
                if not os.path.exists(p):
                    raise colortext.Exception('The folder %s was expected to exist after the preminimization step but could not be found.' % p)

            # Read in the dataset file
            try:
                dataset_filepath = os.path.join(benchmark_run_directory, 'dataset.json')
                dataset = json.loads(read_file(dataset_filepath))
                dataset_cases = {}
                for dsc in dataset['data']:
                    dataset_cases[dsc['RecordID']] = dsc
            except Exception, e:
                raise colortext.Exception('An error occurred parsing the JSON file: %s..' % str(e))

            # Read the previously extracted benchmark data and structural scores (benchmark_data.json) from file or create that file if it does not exist
            benchmark_data_filepath = os.path.join(benchmark_run_directory, 'benchmark_data.json')
            analysis_data = {}
            if not(os.path.exists(benchmark_data_filepath)) or arguments.get('--use_existing_benchmark_data') == 0:
                colortext.warning('Creating %s which contains component and summary scores for each case and generated structure.' % benchmark_data_filepath)
                job_dirs = sorted([jd for jd in glob.glob(os.path.join(ddg_data_dir, '*')) if os.path.isdir(jd) and os.path.split(jd)[1].isdigit()])
                c, num_dirs = 0.0, float(len(job_dirs)) / 100.0
                for jd in job_dirs:
                    jdirname = os.path.split(jd)[1]
                    record_id = int(jdirname)
                    c += 1.0
                    colortext.wcyan('\rProgress: {0:d}%'.format(int(c/num_dirs))); sys.stdout.flush()
                    analysis_data[record_id] = BenchmarkManager.extract_data(jd, self.take_lowest)
                colortext.wlightpurple('\rWriting to file...'); sys.stdout.flush()
                write_file(benchmark_data_filepath, json.dumps(analysis_data, indent = 4, sort_keys=True))
                colortext.wgrey('\r'); sys.stdout.flush()
            else:
                colortext.warning('Found an existing benchmark_data.json file containing component and summary scores for each case and generated structure:')
                colortext.wyellow('\r...loading'); sys.stdout.flush()
                analysis_data_ = json.loads(read_file(benchmark_data_filepath))
                for k, v in analysis_data_.iteritems():
                    analysis_data[int(k)] = v
                colortext.wgrey('\r'); sys.stdout.flush()

            self.benchmark_run_data[benchmark_run_name] = BenchmarkRun(
                benchmark_run_name,
                benchmark_run_directory,
                self.analysis_directory,
                dataset_cases,
                analysis_data,
                arguments['--use_single_reported_value'],
                take_lowest = self.take_lowest,
                generate_plots = self.generate_plots
            )


    def analyze(self):
        '''Runs the analysis for the different benchmark runs.'''

        # change to create analysis directory
        if not(os.path.exists(self.analysis_directory)):
            try:
                os.mkdir(self.analysis_directory)
                assert(os.path.exists(self.analysis_directory))
            except Exception, e:
                raise colortext.Exception('An exception occurred creating the directory %s.' % self.analysis_directory)

        benchmark_runs = sorted(self.benchmark_run_data.keys())
        # todo: add published data here

        for benchmark_run_name, br in sorted(self.benchmark_run_data.iteritems()):
            colortext.message('\nRunning the analysis for {0}.'.format(benchmark_run_name))
            br.create_dataframe()

            pass
            # run individual analysis

        for x in range(len(benchmark_runs)):
            for y in range(x + 1, len(benchmark_runs)):
                pass
                # run comparative analysis


    ###
    #  Data-extraction methods
    ###

    @staticmethod
    def extract_data(ddg_output_path, take_lowest):
        '''Extract the data for one dataset case.
           Note: This function is written to handle data from Rosetta at the time of writing. This will need to be altered to
           work with certain older revisions and may need to be adapted to work with future revisions.'''
        scores = {}
        rosetta_output_lines = BenchmarkManager.read_stdout(ddg_output_path)
        scores['Errors'] = BenchmarkManager.extract_errors(rosetta_output_lines)
        scores['DDG'] = BenchmarkManager.extract_predicted_ddg(rosetta_output_lines)
        scores['DDG_components'] = BenchmarkManager.extract_summary_data(ddg_output_path)
        for k, v in BenchmarkManager.get_ddg_monomer_scores_per_structure(rosetta_output_lines).iteritems():
            assert(k not in scores)
            scores[k] = v

        top_values = set([3, take_lowest])
        for x in top_values:
            scores['DDG_Top%d' % x] = None
            # If there are less mutant scores than wildtype scores or vice versa then only consider as many as the lowest number.
            # This case should rarely occur unless take_lowest is set too high or a small number of structures was created
            fair_top = min(x, len(scores['Mutant_scores']), len(scores['WildType_scores']))
            if fair_top > 0:
                avg_top_mutant_scores = sum(sorted([components['total'] for id, components in scores['Mutant_scores'].iteritems()])[:fair_top])/float(fair_top)
                avg_top_wildtype_scores = sum(sorted([components['total'] for id, components in scores['WildType_scores'].iteritems()])[:fair_top])/float(fair_top)
                scores['DDG_Top%d' % x] = avg_top_mutant_scores - avg_top_wildtype_scores

        return scores


    @staticmethod
    def read_stdout(ddg_output_path):
        try:
            output_file = os.path.join(ddg_output_path, 'rosetta.out.gz')
            rosetta_output = read_file_lines(output_file)
        except:
            raise colortext.Exception('An error occurred reading the output file %s.' % output_file)
        return rosetta_output


    @staticmethod
    def extract_errors(rosetta_output_lines):
        '''Reads in errors from the predictions which can be useful to report to the user.'''
        errors = {
            'Derivative error count' : 0,
        }
        derivative_error_count = 0
        for l in rosetta_output_lines:
            if l.lower().find('inaccurate g') != -1:
                derivative_error_count += 1
        errors['Derivative error count'] = derivative_error_count
        return errors


    @staticmethod
    def extract_predicted_ddg(rosetta_output_lines):
        ddGMover_lines = [l for l in rosetta_output_lines if l.strip() and l.startswith("protocols.moves.ddGMover: mutate")]
        assert(len(ddGMover_lines) == 1)
        ddGMover_line = ddGMover_lines[0].strip()
        mtchs = ddGMover_regex.match(ddGMover_line)
        assert(mtchs)
        return float(mtchs.group(1))


    @staticmethod
    def extract_summary_data(ddg_output_path):
        score_data = [l for l in read_file_lines(os.path.join(ddg_output_path, 'ddg_predictions.out')) if l.strip()]
        assert(len(score_data) == 2) # Assuming only one line here

        score_headers = score_data[0].split()
        score_data = score_data[1].split()

        score_summary = {}
        assert(len(score_data) > 2 and len(score_data) == len(score_headers))
        assert(score_data[0] == "ddG:" and score_headers[0] == "ddG:")
        for x in range(1, len(score_headers)):
            try:
                score_summary[score_headers[x]] = float(score_data[x])
            except:
                score_summary[score_headers[x]] = score_data[x]
        return score_summary


    @staticmethod
    def get_ddg_monomer_scores_per_structure(rosetta_output_lines):
        '''Returns a dict mapping the DDG scores from a ddg_monomer run to a list of structure numbers.'''

        # Parse the stdout into two mappings (one for wildtype structures, one for mutant structures) mapping
        # structure IDs to a dict containing the score components
        wildtype_scores = {}
        mutant_scores = {}
        s1 = 'score before mutation: residue'
        s1_len = len(s1)
        s2 = 'score after mutation: residue'
        s2_len = len(s2)
        for line in rosetta_output_lines:
            idx = line.find(s1)
            if idx != -1:
                idx += s1_len
                mtchs = re.match('.*?(\d+) %s' % s1, line)
                structure_id = int(mtchs.group(1))
                assert(structure_id not in wildtype_scores)
                tokens = line[idx:].split()
                d = {'total' : float(tokens[0])}
                for x in range(1, len(tokens), 2):
                    component_name = tokens[x].replace(':', '')
                    component_value = float(tokens[x + 1])
                    d[component_name] = component_value
                wildtype_scores[structure_id] = d
            else:
                idx = line.find(s2)
                if idx != -1:
                    idx += s2_len
                    mtchs = re.match('.*?(\d+) %s' % s2, line)
                    structure_id = int(mtchs.group(1))
                    assert(structure_id not in mutant_scores)
                    tokens = line[idx:].split()
                    d = {'total' : float(tokens[1])}
                    for x in range(2, len(tokens), 2):
                        component_name = tokens[x].replace(':', '')
                        component_value = float(tokens[x + 1])
                        d[component_name] = component_value
                    mutant_scores[structure_id] = d

        # Sanity checks
        num_structures = max(wildtype_scores.keys())
        expected_keys = set(range(1, num_structures + 1))
        assert(expected_keys == set(wildtype_scores.keys()))
        assert(expected_keys == set(mutant_scores.keys()))

        # Create a list of lists - MutantScoreOrder - of structure IDs e.g. [[5,1,34], [23], [12,3], ...] which is ordered
        # by increasing energy so that each sublist contains structure IDs of equal energy and if structures have the same
        # energy then their IDs are in the same sublist
        d = {}
        for structure_id, scores in sorted(mutant_scores.iteritems()):
            d[scores['total']] = d.get(scores['total'], [])
            d[scores['total']].append(structure_id)
        MutantScoreOrder = []
        for score, structure_ids in sorted(d.iteritems()):
            MutantScoreOrder.append(structure_ids)

        # Sanity check - make sure that MutantScoreOrder is really ordered such that each set of structure IDs contains
        # structures of the same energy and of a lower energy than the following set of structure IDs in the list
        for x in range(len(MutantScoreOrder) - 1):
            s1 = set([mutant_scores[n]['total'] for n in MutantScoreOrder[x]])
            assert(len(s1) == 1)
            if x + 1 < len(MutantScoreOrder):
                s2 = set([mutant_scores[n]['total'] for n in MutantScoreOrder[x + 1]])
                assert(len(s2) == 1)
                assert(s1.pop() < s2.pop())

        # Determine the lowest scoring pair of structures
        # Iterate over the (wildtype, mutant) pairs and determine the structure ID for the pair with the lowest energy mutant
        # Note: There are multiple ways to select the best pair. For example, if multiple mutants have the same minimal total
        # score, we could have multiple wildtype structures to choose from. In this case, we choose a pair where the wildtype
        # structure has the minimal total score.

        lowest_mutant_score = min([v['total'] for k, v in mutant_scores.iteritems()])
        mutant_structure_ids  = [k for k, v in mutant_scores.iteritems() if v['total'] == lowest_mutant_score]
        if len(mutant_structure_ids) > 1:
            lowest_wildtype_score = min([v['total'] for k, v in wildtype_scores.iteritems() if k in mutant_structure_ids])
            lowest_scoring_pair_id = sorted([k for k, v in wildtype_scores.iteritems() if v['total'] == lowest_wildtype_score])[0] # this is deterministic but there could potentially be multiple valid IDs
        else:
            lowest_scoring_pair_id = mutant_structure_ids[0]

        return dict(
            Lowest_scoring_pair = lowest_scoring_pair_id,
            WildType_scores = wildtype_scores,
            Mutant_scores = mutant_scores,
            Mutant_score_order = MutantScoreOrder,
        )



class BenchmarkRun(object):
    '''A object to contain benchmark run data which can be used to analyze that run or else to cross-analyze the run with another run.'''

    # Class variables
    amino_acid_details = {}
    CAA, PAA, HAA = set(), set(), set()


    def __init__(self, benchmark_run_name, benchmark_run_directory, analysis_directory, dataset_cases, analysis_data, use_single_reported_value, take_lowest = 3, generate_plots = True):
        self.amino_acid_details, self.CAA, self.PAA, self.HAA = BenchmarkRun.get_amino_acid_details()
        self.benchmark_run_name = benchmark_run_name
        self.benchmark_run_directory = benchmark_run_directory
        self.dataset_cases = dataset_cases
        self.analysis_data = analysis_data
        self.analysis_directory = analysis_directory
        self.analysis_file_prefix = os.path.join(self.analysis_directory, self.benchmark_run_name + '_')
        print(self.analysis_file_prefix)
        self.use_single_reported_value = use_single_reported_value
        self.generate_plots = generate_plots
        self.take_lowest = take_lowest


    @staticmethod
    def get_amino_acid_details():
        if not BenchmarkRun.amino_acid_details:
            # Amino acid properties
            polarity_map = {'polar' : 'P', 'charged' : 'C', 'hydrophobic' : 'H'}
            aromaticity_map = {'aliphatic' : 'L', 'aromatic' : 'R', 'neither' : '-'}
            amino_acid_detail_headers = 'Code,Long code,Name,Polarity,Aromaticity,Hydrophobicity pH7,Sidechain acidity,pKa,Average mass,van der Waals volume,Size,Is tiny?'
            amino_acid_details_ = [
                'A,ALA,Alanine,non-polar,aliphatic,hydrophobic,neutral,NULL,71.0788,67,small,1',
                'C,CYS,Cysteine,polar,neither,hydrophilic,neutral,8.7,103.1388,86,small,1',
                'D,ASP,Aspartic acid,charged,neither,hydrophilic,acidic,3.9,115.0886,91,small,0',
                'E,GLU,Glutamic acid,charged,neither,hydrophilic,acidic,4.5,129.1155,109,large,0',
                'F,PHE,Phenylalanine,non-polar,aromatic,hydrophobic,neutral,NULL,147.1766,135,large,0',
                'G,GLY,Glycine,polar,neither,hydrophilic,neutral,NULL,57.0519,48,small,1',
                'H,HIS,Histidine,charged,neither,hydrophilic,basic,6.04,137.1411,118,large,0',
                'I,ILE,Isoleucine,non-polar,aliphatic,hydrophobic,neutral,NULL,113.1594,124,large,0',
                'K,LYS,Lysine,charged,neither,hydrophilic,basic,10.54,128.1741,135,large,0',
                'L,LEU,Leucine,non-polar,aliphatic,hydrophobic,neutral,NULL,113.1594,124,large,0',
                'M,MET,Methionine,non-polar,aliphatic,hydrophobic,neutral,NULL,131.1986,124,large,0',
                'N,ASN,Asparagine,polar,neither,hydrophilic,neutral,NULL,114.1039,96,small,0',
                'P,PRO,Proline,non-polar,neither,hydrophobic,neutral,NULL,97.1167,90,small,0',
                'Q,GLN,Glutamine,polar,neither,hydrophilic,neutral,NULL,128.1307,114,large,0',
                'R,ARG,Arginine,charged,neither,hydrophilic,basic,12.48,156.1875,148,large,0',
                'S,SER,Serine,polar,neither,hydrophilic,neutral,NULL,87.0782,73,small,1',
                'T,THR,Threonine,polar,neither,hydrophilic,neutral,NULL,101.1051,93,small,0',
                'V,VAL,Valine,non-polar,aliphatic,hydrophobic,neutral,NULL,99.1326,105,small,0',
                'W,TRP,Tryptophan,non-polar,aromatic,hydrophobic,neutral,NULL,186.2132,163,large,0',
                'Y,TYR,Tyrosine,polar,aromatic,hydrophobic,neutral,10.46,163.176,141,large,0' # Note: we treat tyrosine as hydrophobic in the polar/charged vs hydrophobic/Non-polar plot
            ]

            amino_acid_detail_headers = [t.strip() for t in amino_acid_detail_headers.split(',') if t.strip()]
            for aad in amino_acid_details_:
                tokens = aad.split(',')
                assert(len(tokens) == len(amino_acid_detail_headers))
                d = {}
                for x in range(len(amino_acid_detail_headers)):
                    d[amino_acid_detail_headers[x]] = tokens[x]
                aa_code = d['Code']
                BenchmarkRun.amino_acid_details[aa_code] = d
                del d['Code']
                d['Polarity'] = polarity_map.get(d['Polarity'], 'H')
                d['Aromaticity'] = aromaticity_map[d['Aromaticity']]
                d['Average mass'] = float(d['Average mass'])
                d['Is tiny?'] = d['Is tiny?'] == 1
                d['van der Waals volume'] = float(d['van der Waals volume'])
                try: d['pKa'] = float(d['pKa'])
                except: d['pKa'] = None

                if aa_code == 'Y':
                    BenchmarkRun.HAA.add(aa_code) # Note: Treating tyrosine as hydrophobic
                elif d['Polarity'] == 'C':
                    BenchmarkRun.CAA.add(aa_code)
                elif d['Polarity'] == 'P':
                    BenchmarkRun.PAA.add(aa_code)
                elif d['Polarity'] == 'H':
                    BenchmarkRun.HAA.add(aa_code)
        assert(len(BenchmarkRun.CAA.intersection(BenchmarkRun.PAA)) == 0 and len(BenchmarkRun.PAA.intersection(BenchmarkRun.HAA)) == 0 and len(BenchmarkRun.HAA.intersection(BenchmarkRun.CAA)) == 0)
        return BenchmarkRun.amino_acid_details, BenchmarkRun.CAA, BenchmarkRun.PAA, BenchmarkRun.HAA


    def create_dataframe(self):

        # Create XY data
        analysis_csv_input_filepath = os.path.join(self.benchmark_run_directory, 'analysis_input.csv')
        analysis_json_input_filepath = os.path.join(self.benchmark_run_directory, 'analysis_input.json')
        analysis_pandas_input_filepath = os.path.join(self.benchmark_run_directory, 'analysis_input.pandas')
        analysis_data = self.analysis_data
        dataset_cases = self.dataset_cases

        print('Creating the analysis input file %s and human-readable CSV and JSON versions %s and %s.' % (analysis_pandas_input_filepath, analysis_csv_input_filepath, analysis_json_input_filepath))
        if len(analysis_data) > len(dataset_cases):
            raise colortext.Exception('ERROR: There seems to be an error - there are more predictions than cases in the dataset. Exiting.')
        elif len(analysis_data) < len(dataset_cases):
            colortext.error('\nWARNING: %d cases missing for analysis; there are %d predictions in the output directory but %d cases in the dataset. The analysis below does not cover the complete dataset.\n' % (len(dataset_cases) - len(analysis_data), len(analysis_data), len(dataset_cases)))

        # ddg_analysis_type can be set to 'DDG' or 'DDG_Top[x]' (e.g. 'DDG_Top3').
        # 'DDG' uses the value reported by ddg_monomer.
        # 'DDG_Top3' (generated by default) uses the metric from Kellogg et al. based on the three lowest scoring mutant structures and the three lowest scoring wildtype structures
        ddg_analysis_type = 'DDG_Top%d' % self.take_lowest
        if arguments['--use_single_reported_value']:
            ddg_analysis_type = 'DDG'
            print('\nThe predicted DDG value per case is the single DDG value reported by the ddg_monomer application.')
        elif self.take_lowest == 3:
            print('\nThe predicted DDG value per case is computed using the {0} lowest-scoring mutant structures and the {0} lowest-scoring wildtype structures as in the paper by Kellogg et al.'.format(self.take_lowest))
        else:
            print('\nThe predicted DDG value per case is computed using the {0} lowest-scoring mutant structures and the {0} lowest-scoring wildtype structures.'.format(self.take_lowest))
        sys.exit(0)
        json_records = dict(
            All = [],
                ByVolume = dict(
                    SL = [], LS = [], XX = [], O = []
                ),
                GP = {
                    True : [], False : []
                },
        )
        csv_headers=[
            'DatasetID',
            'PDBFileID',
            'Mutations',
            'Experimental',
            'Predicted',
            'AbsoluteError',
            'StabilityClassification',
            'ResidueCharges',
            'VolumeChange',
            'WildTypeDSSPType',
            'WildTypeDSSPSimpleSSType',
            'WildTypeDSSPExposure',
            'WildTypeSCOPClass',
            'WildTypeSCOPFold',
            'WildTypeSCOPClassification',
            'WildTypeExposure',
            'WildTypeAA',
            'MutantAA',
            'PDBResolution',
            'PDBResolutionBin',
            'MonomerLength',
        ]
        csv_file = [','.join(csv_headers)]

        SCOP_classifications = set()
        SCOP_folds = set()
        SCOP_classes = set()

        include_derived_mutations = arguments['--include_derived_mutations']
        for record_id, predicted_data in sorted(analysis_data.iteritems()):
            record = dataset_cases[record_id]
            if record['DerivedMutation'] and not include_derived_mutations:
                continue

            full_scop_classification, scop_class, scop_fold = None, None, None
            mutations = record['Mutations']
            scops = set([m['SCOP class'] for m in mutations])
            if len(scops) > 1:
                colortext.warning('Warning: There is more than one SCOPe class for record {0}.'.format(record_id))
            else:
                full_scop_classification = scops.pop()
                scop_tokens = full_scop_classification.split('.')
                scop_class = scop_tokens[0]
                if len(scop_tokens) > 1:
                    scop_fold = '.'.join(scop_tokens[0:2])

                SCOP_classifications.add(full_scop_classification)
                SCOP_classes.add(scop_class)
                SCOP_folds.add(scop_fold)

            DSSPSimpleSSType, DSSPSimpleSSTypes = None, set([m['DSSPSimpleSSType'] for m in mutations])
            if len(DSSPSimpleSSTypes) == 1: DSSPSimpleSSType = DSSPSimpleSSTypes.pop()

            DSSPType, DSSPTypes = None, set([m['DSSPType'] for m in mutations])
            if len(DSSPTypes) == 1: DSSPType = DSSPTypes.pop()

            DSSPExposure, DSSPExposures = None, set([m['DSSPExposure'] for m in mutations])
            if len(DSSPExposures) == 1: DSSPExposure = DSSPExposures.pop()

            mutation_string = '; '.join(['{0} {1}{2}{3}'.format(m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA']) for m in mutations])

            residue_charge, residue_charges = None, set()
            exposure, exposures = None, set()
            volume_change, volume_changes = None, set()
            for m in mutations:

                wtaa = m['WildTypeAA']
                mutaa = m['MutantAA']

                # Burial
                if m['DSSPExposure'] != None:
                    if m['DSSPExposure'] > burial_cutoff:
                        exposures.add('E')
                    else:
                        exposures.add('B')
                else:
                    exposures.add(None)

                # Volume
                if amino_acid_details[wtaa]['van der Waals volume'] < amino_acid_details[mutaa]['van der Waals volume']:
                    volume_changes.add('SL')
                elif amino_acid_details[wtaa]['van der Waals volume'] > amino_acid_details[mutaa]['van der Waals volume']:
                    volume_changes.add('LS')
                elif amino_acid_details[wtaa]['van der Waals volume'] == amino_acid_details[mutaa]['van der Waals volume']:
                    volume_changes.add('XX')

                # Charge
                if ((wtaa in CAA or wtaa in PAA) and (mutaa in HAA)) or ((mutaa in CAA or mutaa in PAA) and (wtaa in HAA)):
                    residue_charges.add('Change')
                elif (wtaa in CAA or wtaa in PAA) and (mutaa in CAA or mutaa in PAA):
                    residue_charges.add('Polar/Charged')
                elif (wtaa in HAA) and (mutaa in HAA):
                    residue_charges.add('Hydrophobic/Non-polar')
                else:
                     raise colortext.Exception('Should not reach here.')

            record_wtaa, wtaas = None, set([m['WildTypeAA'] for m in mutations])
            if len(wtaas) == 1: record_wtaa = wtaas.pop()

            record_mutaa, mutaas = None, set([m['MutantAA'] for m in mutations])
            if len(mutaas) == 1: record_mutaa = mutaas.pop()

            pdb_chain, pdb_chains = None, set([m['Chain'] for m in mutations])
            assert(len(pdb_chains) == 1)
            pdb_chain = pdb_chains.pop()

            # Take unique values
            if len(residue_charges) == 1: residue_charge = residue_charges.pop()
            if len(exposures) == 1: exposure = exposures.pop()
            if len(volume_changes) == 1: volume_change = volume_changes.pop()

            # Set the PDB input path
            pdb_data = {}
            try:
                pdb_data_ = json.loads(read_file('../../input/json/pdbs.json'))
                for k, v in pdb_data_.iteritems():
                    pdb_data[k.upper()] = v
            except Exception, e:
                colortext.error('input/json/pdbs.json could not be found - PDB-specific analysis cannot be performed.')

            stability_classification = fraction_correct([record['DDG']], [predicted_data[ddg_analysis_type]], x_cutoff = stability_classication_x_cutoff, y_cutoff = stability_classication_y_cutoff)
            assert(stability_classification == 0 or stability_classification == 1)
            stability_classification = int(stability_classification)

            # Partition the data by PDB resolution with bins: N/A, <1.5, 1.5-<2.0, 2.0-<2.5, >=2.5
            pdb_record = pdb_data.get(record['PDBFileID'].upper())
            pdb_resolution_bin = None
            pdb_resolution = pdb_record.get('Resolution')
            if pdb_resolution != None:
                if pdb_resolution < 1.5:
                    pdb_resolution_bin = '<1.5'
                elif pdb_resolution < 2.0:
                    pdb_resolution_bin = '1.5-2.0'
                elif pdb_resolution < 2.5:
                    pdb_resolution_bin = '2.0-2.5'
                else:
                    pdb_resolution_bin = '>=2.5'
            pdb_resolution_bin = pdb_resolution_bin or 'N/A'

            # Create the data matrix
            json_record = dict(
                DatasetID = record_id,
                PDBFileID = record['PDBFileID'],
                Mutations = mutation_string,
                Experimental = record['DDG'],
                Predicted = predicted_data[ddg_analysis_type],
                AbsoluteError = abs(record['DDG'] - predicted_data[ddg_analysis_type]),
                StabilityClassification = stability_classification,
                ResidueCharges = residue_charge,
                VolumeChange = volume_change,
                WildTypeDSSPType = DSSPType,
                WildTypeDSSPSimpleSSType = DSSPSimpleSSType,
                WildTypeDSSPExposure = DSSPExposure,
                WildTypeSCOPClass = scop_class,
                WildTypeSCOPFold = scop_fold,
                WildTypeSCOPClassification = full_scop_classification,
                WildTypeExposure = exposure,
                WildTypeAA = record_wtaa,
                MutantAA = record_mutaa,
                PDBResolution = pdb_record.get('Resolution'),
                PDBResolutionBin = pdb_resolution_bin,
                MonomerLength = len(pdb_record.get('Chains', {}).get(pdb_chain, {}).get('Sequence', '')) or None,
                )

            json_records['ByVolume'][determine_SL_class(record)].append(json_record)
            json_records['GP'][has_G_or_P(record)].append(json_record)
            json_records['All'].append(json_record)

            for h in csv_headers:
                assert(',' not in str(json_record[h]))
            csv_file.append(','.join([str(json_record[h]) for h in csv_headers]))

        colortext.message('The mutated residues span {0} unique SCOP(e) classifications in {1} unique SCOP(e) folds and {2} unique SCOP(e) classes.'.format(len(SCOP_classifications), len(SCOP_folds), len(SCOP_classes)))

        write_file(analysis_csv_input_filepath, '\n'.join(csv_file))
        write_file(analysis_json_input_filepath, json.dumps(json_records['All'], indent = 4, sort_keys=True))

        ##todo self.generate_plots

    def compare(self, other):
        '''Compare this benchmark run with another run.'''
        pass
        #todo implement this


    def analyze_run(self):
        '''This function runs the analysis, creating plots and the summary file.'''

        rel_plot_filename_prefix = arguments['--plot_filename_prefix'][0]
        subplot_dir = os.path.join(output_dir, 'subplots')
        try: os.mkdir(subplot_dir)
        except: pass
        if not os.path.exists(subplot_dir):
            raise colortext.Exception('Could not create subdirectory {0} for the plots.'.format(subplot_dir))
        plot_filename_prefix = os.path.join(subplot_dir, '{0}'.format(rel_plot_filename_prefix))

        # Plot the optimum y-cutoff over a range of x-cutoffs for the fraction correct metric. Include the user's cutoff in the range
        scalar_adjustment = plot_optimum_prediction_fraction_correct_cutoffs_over_range(plot_filename_prefix, json_records['All'], min(stability_classication_x_cutoff, 0.5), max(stability_classication_x_cutoff, 3.0))

        # Plot which y-cutoff yields the best value for the fraction correct metric
        plot_optimum_prediction_fraction_correct_cutoffs(plot_filename_prefix, json_records['All'], stability_classication_x_cutoff)

        # Create an adjusted set of records scaled according to the fraction correct fitting
        adjusted_records = copy.deepcopy(json_records['All'])
        for record in adjusted_records:
            record['Predicted'] = record['Predicted'] / scalar_adjustment

        # Create a scatterplot for the results
        output_filename = '{0}_scatterplot.png'.format(plot_filename_prefix)
        print('Saving scatterplot to %s.' % output_filename)
        plot(json_records['All'], output_filename, RInterface.correlation_coefficient_gplot, title = 'Experimental vs. Prediction')

        # Create a scatterplot for the adjusted results
        output_filename = '{0}_scatterplot_adjusted_with_scalar.png'.format(plot_filename_prefix)
        print('Saving scatterplot to %s.' % output_filename)
        plot(adjusted_records, output_filename, RInterface.correlation_coefficient_gplot, title = 'Experimental vs. Prediction: adjusted scale')

        # Plot a histogram of the absolute errors
        plot_absolute_error_histogram(plot_filename_prefix, json_records['All'])

        # Create a series of scatterplots colored by different criteria
        scatterplot_generic('Experimental vs. Prediction - Exposure (cutoff = %0.2f)' % burial_cutoff, json_records['All'], scatterplot_exposure, '{0}_scatterplot_exposure.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Residue charges', json_records['All'], scatterplot_charges, '{0}_scatterplot_charges.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Change in volume', json_records['All'], scatterplot_volume, '{0}_scatterplot_volume.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Wildtype residue s.s.', json_records['All'], scatterplot_ss, '{0}_scatterplot_ss.png'.format(plot_filename_prefix))
        if len(SCOP_classes) <= 25:
            scatterplot_generic('Experimental vs. Prediction - WT residue SCOP class', json_records['All'], scatterplot_scop_class, '{0}_scatterplot_scop_class.png'.format(plot_filename_prefix))
        if len(SCOP_folds) <= 25:
            scatterplot_generic('Experimental vs. Prediction - WT residue SCOP fold', json_records['All'], scatterplot_scop_fold, '{0}_scatterplot_scop_fold.png'.format(plot_filename_prefix))
        if len(SCOP_classifications) <= 25:
            scatterplot_generic('Experimental vs. Prediction - WT residue SCOP classification', json_records['All'], scatterplot_scop_classification, '{0}_scatterplot_scop_classification.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Wildtype', json_records['All'], scatterplot_wildtype_aa, '{0}_scatterplot_wildtype_aa.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Mutant', json_records['All'], scatterplot_mutant_aa, '{0}_scatterplot_mutant_aa.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Glycine/Proline', json_records['All'], scatterplot_GP, '{0}_scatterplot_gp.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - PDB resolution', json_records['All'], scatterplot_pdb_res_binned, '{0}_scatterplot_pdb_res_binned.png'.format(plot_filename_prefix))
        scatterplot_generic('Experimental vs. Prediction - Chain length', json_records['All'], scatterplot_chain_length, '{0}_scatterplot_chain_length.png'.format(plot_filename_prefix))

        if include_derived_mutations:
            colortext.message('\nRunning analysis (derived mutations will be included):')
        else:
            colortext.message('\nRunning analysis (derived mutations will be omitted):')
        colortext.warning('The stability classification cutoffs are: Experimental %0.2f kcal/mol, Predicted: %0.2f energy units.' % (stability_classication_x_cutoff, stability_classication_y_cutoff))

        volume_groups = {}
        for aa_code, aa_details in amino_acid_details.iteritems():
            v = int(aa_details['van der Waals volume']) # Note: I only convert to int here to match the old script behavior and because all volumes are integer values so it does not do any harm
            volume_groups[v] = volume_groups.get(v, [])
            volume_groups[v].append(aa_code)

        print('\n\nSection 1. Breakdown by volume.')
        print('A case is considered a small-to-large (resp. large-to-small) mutation if all of the wildtype residues have a smaller (resp. larger) van der Waals volume than the corresponding mutant residue. The order is defined as %s so some cases are considered to have no change in volume e.g. MET -> LEU.' % (' < '.join([''.join(sorted(v)) for k, v in sorted(volume_groups.iteritems())])))
        for subcase in ('XX', 'SL', 'LS'):
            print('\n' + '*'*10 + (' Statistics - %s (%d cases)' % (by_volume_descriptions[subcase], len(json_records['ByVolume'][subcase]))) +'*'*10)
            print(format_stats_for_printing(get_xy_dataset_statistics(json_records['ByVolume'][subcase], fcorrect_x_cutoff = stability_classication_x_cutoff, fcorrect_y_cutoff = stability_classication_y_cutoff)))

        print('\n\nSection 2. Separating out mutations involving glycine or proline.')
        print('This cases may involve changes to secondary structure so we separate them out here.')
        print('\n' + '*'*10 + (' Statistics - cases with G or P (%d cases)' % (len(json_records['GP'][True]))) +'*'*10)
        print(format_stats_for_printing(get_xy_dataset_statistics(json_records['GP'][True], fcorrect_x_cutoff = stability_classication_x_cutoff, fcorrect_y_cutoff = stability_classication_y_cutoff)))
        print('\n' + '*'*10 + (' Statistics - cases without G or P (%d cases)' % (len(json_records['GP'][False]))) +'*'*10)
        print(format_stats_for_printing(get_xy_dataset_statistics(json_records['GP'][False], fcorrect_x_cutoff = stability_classication_x_cutoff, fcorrect_y_cutoff = stability_classication_y_cutoff)))

        print('\nSection 3. Entire dataset using a scaling factor of 1/%.03f to improve the fraction correct metric.' % scalar_adjustment)
        colortext.warning('Warning: Results in this section use an averaged scaling factor to improve the value for the fraction correct metric. This scalar will vary over benchmark runs so these results should not be interpreted as performance results; they should be considered as what could be obtained if the predicted values were scaled by a "magic" value.')
        print('\n' + '*'*10 + (' Statistics - complete dataset (%d cases)' % len(adjusted_records)) +'*'*10)
        # For these statistics, we assume that we have reduced any scaling issues and use the same cutoff for the Y-axis as the user specified for the X-axis
        print(format_stats_for_printing(get_xy_dataset_statistics(adjusted_records, fcorrect_x_cutoff = stability_classication_x_cutoff, fcorrect_y_cutoff = stability_classication_x_cutoff)))

        print('\n\nSection 4. Entire dataset.')
        colortext.message('Overall statistics')
        print('\n' + '*'*10 + (' Statistics - complete dataset (%d cases)' % len(json_records['All'])) +'*'*10)
        print(format_stats_for_printing(get_xy_dataset_statistics(json_records['All'], fcorrect_x_cutoff = stability_classication_x_cutoff, fcorrect_y_cutoff = stability_classication_y_cutoff)))

        # Combine the plots into a PDF file
        plot_file = os.path.join(output_dir, rel_plot_filename_prefix + '_benchmark_plots.pdf')
        colortext.message('\n\nCreating a PDF containing all of the plots: {0}'.format(plot_file))
        try:
            if os.path.exists(plot_file):
                os.remove(plot_file)
            p = subprocess.Popen(shlex.split('convert {0} {1}'.format(os.path.join(subplot_dir, '*.png'), rel_plot_filename_prefix + '_benchmark_plots.pdf')), cwd = output_dir)
            #print(shlex.split('convert {0} {1}'.format(os.path.join(subplot_dir, '*.png'), rel_plot_filename_prefix + '_benchmark_plots.pdf')))
            stdoutdata, stderrdata = p.communicate()
            if p.returncode != 0: raise Exception('')
        except:
            colortext.error('An error occurred while combining the positional scatterplots using the convert application (ImageMagick).')


    def determine_optimum_fraction_correct_cutoffs(self, records, stability_classication_x_cutoff):
        '''Determines the value of stability_classication_y_cutoff which approximately maximizes the fraction correct
           measurement w.r.t. a fixed stability_classication_x_cutoff. This function uses discrete sampling and so it
           may miss the actual maximum. We use two rounds of sampling: i) a coarse-grained sampling (0.1 energy unit
           intervals); and ii) finer sampling (0.01 unit intervals).
           In both rounds, we choose the one corresponding to a lower value for the cutoff in cases of multiple maxima.'''

        # Determine the value for the fraction correct y-value (predicted) cutoff which will approximately yield the
        # maximum fraction-correct value
        x_values = [record['Experimental'] for record in records]
        y_values = [record['Predicted'] for record in records]
        fraction_correct_range = []

        # Round 1 : Coarse sampling. Test 0.5 -> 8.0 in 0.1 increments
        for z in range(5, 80):
            w = float(z) / 10.0
            fraction_correct_range.append((w, fraction_correct(x_values, y_values, x_cutoff = stability_classication_x_cutoff, y_cutoff = w)))
        max_value_cutoff, max_value = fraction_correct_range[0][0], fraction_correct_range[0][1]
        for p in fraction_correct_range:
            if p[1] > max_value:
                max_value_cutoff, max_value = p[0], p[1]

        # Round 2 : Finer sampling. Test max_value_cutoff - 0.1 -> max_value_cutoff + 0.1 in 0.01 increments
        for z in range(int((max_value_cutoff - 0.1) * 100), int((max_value_cutoff + 0.1) * 100)):
            w = float(z) / 100.0
            fraction_correct_range.append((w, fraction_correct(x_values, y_values, x_cutoff = stability_classication_x_cutoff, y_cutoff = w)))
        fraction_correct_range = sorted(set(fraction_correct_range)) # sort so that we find the lowest cutoff value in case of duplicate fraction correct values
        max_value_cutoff, max_value = fraction_correct_range[0][0], fraction_correct_range[0][1]
        for p in fraction_correct_range:
            if p[1] > max_value:
                max_value_cutoff, max_value = p[0], p[1]

        return max_value_cutoff, max_value, fraction_correct_range


    def plot_optimum_prediction_fraction_correct_cutoffs(self, plot_filename_prefix, records, stability_classication_x_cutoff):

        # Determine the optimal values
        max_value_cutoff, max_value, fraction_correct_range = determine_optimum_fraction_correct_cutoffs(records, stability_classication_x_cutoff)

        # Filenames
        output_filename_prefix = '{0}_optimum_fraction_correct_at_{1}_kcal_mol'.format(plot_filename_prefix, '%.2f' % stability_classication_x_cutoff)
        plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'
        print('Saving scatterplot to %s.' % plot_filename)

        # Create CSV input
        lines = ['NeutralityCutoff,FractionCorrect,C']
        for p in fraction_correct_range:
            if p[1] == max_value:
                lines.append(','.join(map(str, (p[0], p[1], 'best'))))
            else:
                lines.append(','.join(map(str, (p[0], p[1], 'other'))))
        write_file(csv_filename, '\n'.join(lines))

        # Create plot
        title = 'Optimum cutoff for fraction correct metric at %0.2f kcal/mol' % stability_classication_x_cutoff
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

plot_scale <- scale_color_manual(
values = c( "best" = '#00dd00', "other" = '#666666'),
labels = c( "best" = "Best", "other" = "Other"),
guide = "none") # do not show the legend

best_y = max(plot_data$FractionCorrect)
p <- ggplot(data = plot_data, aes(x = NeutralityCutoff, y = FractionCorrect)) +
 plot_scale +
 xlab("Neutrality cutoff (energy units)") +
 ylab("Fraction correct") +
 ggtitle("%(title)s") +
 geom_point(aes(color = C)) +
 geom_line() +
 geom_smooth() +
 geom_text(hjust=0, size=4, color="black", aes(6.5, best_y, fontface="plain", family = "sans", label=sprintf("Max = %(max_value)0.2f\\nCutoff = %(max_value_cutoff)0.2f")))
p
dev.off()'''
        RInterface._runRScript(r_script % locals())


    def plot_optimum_prediction_fraction_correct_cutoffs_over_range(self, plot_filename_prefix, records, min_stability_classication_x_cutoff, max_stability_classication_x_cutoff):
        '''Plots the optimum cutoff for the predictions to maximize the fraction correct metric over a range of experimental cutoffs.
           Returns the average scalar corresponding to the best value of fraction correct over a range of cutoff values for the experimental cutoffs.'''

        # Filenames
        output_filename_prefix = '{0}_optimum_fraction_correct_at_varying_kcal_mol'.format(plot_filename_prefix)
        plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'
        print('Saving scatterplot to %s.' % plot_filename)

        # Create CSV input
        lines = ['ExperimentalCutoff,BestPredictionCutoff']
        x_cutoff = min_stability_classication_x_cutoff
        x_values = []
        y_values = []
        avg_scale = 0
        while x_cutoff < max_stability_classication_x_cutoff + 0.1:
            max_value_cutoff, max_value, fraction_correct_range = determine_optimum_fraction_correct_cutoffs(records, x_cutoff)
            lines.append(','.join(map(str, (x_cutoff, max_value_cutoff))))
            x_values.append(x_cutoff)
            y_values.append(max_value_cutoff)
            avg_scale += max_value_cutoff / x_cutoff
            x_cutoff += 0.1
        write_file(csv_filename, '\n'.join(lines))

        # Determine the average scalar needed to fit the plot
        avg_scale = avg_scale / len(x_values)
        x_values = numpy.array(x_values)
        y_values = numpy.array(y_values)
        scalars = y_values / x_values
        average_scalar = numpy.mean(scalars)
        plot_label_1 = 'Scalar == %0.2f' % average_scalar
        plot_label_2 = 'sigma == %0.2f' % numpy.std(scalars)

        # Create plot
        title = 'Optimum cutoff for fraction correct metric at varying experimental cutoffs'
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

max_y = max(plot_data$BestPredictionCutoff)
p <- ggplot(data = plot_data, aes(x = ExperimentalCutoff, y = BestPredictionCutoff)) +
 plot_scale +
 xlab("Experimental cutoff (kcal/mol)") +
 ylab("Optimal prediction cutoff (energy units)") +
 ggtitle("%(title)s") +
 geom_point() +
 geom_line() +
 geom_smooth() +
 geom_text(hjust=0, size=4, color="black", aes(0.5, max_y, fontface="plain", family = "sans", label="%(plot_label_1)s"), parse = T) +
 geom_text(hjust=0, size=4, color="black", aes(0.5, max_y - 0.5, fontface="plain", family = "sans", label="%(plot_label_2)s"), parse = T)
p
dev.off()'''
        RInterface._runRScript(r_script % locals())
        return average_scalar


    def plot_absolute_error_histogram(self, plot_filename_prefix,data):

        # Filenames
        output_filename_prefix = '{0}_absolute_errors'.format(plot_filename_prefix)
        plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'
        print('Saving scatterplot to %s.' % plot_filename)

        # Create CSV input
        lines = ['DatasetID,AbsoluteError']
        for record in data:
            lines.append('{0},{1}'.format(record['DatasetID'], record['AbsoluteError']))
        write_file(csv_filename, '\n'.join(lines))

        # Create plot
        title = 'Distribution of absolute errors (prediction - observed)'
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

m <- ggplot(plot_data, aes(x=AbsoluteError)) +
    geom_histogram(colour = "darkgreen", fill = "green", binwidth = 0.5) +
    ggtitle("%(title)s") +
    xlab("Absolute error (kcal/mol - energy units)") +
    ylab("Number of cases")
m
dev.off()'''
        RInterface._runRScript(r_script % locals())


    def scatterplot_generic(self, title, data, plotfn, plot_filename):
        csv_filename = os.path.splitext(plot_filename)[0] + '.txt'
        plot_commands = plotfn(data, title, csv_filename)
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

%(plot_commands)s

dev.off()''' % locals()
        RInterface._runRScript(r_script)


    def scatterplot_color_by_series(self, data, colorseries, xseries = "Experimental", yseries = "Predicted", title = '', plot_scale = '', point_opacity = 0.4, extra_commands = ''):

        mae = sum([abs(record[xseries] - record[yseries]) for record in data]) / len(data)
        plot_scale_line = ''
        plot_scale_argument = ''
        if plot_scale:
            plot_scale_line = plot_scale.strip()
            plot_scale_argument = '\n    plot_scale +'

        return '''
opacity <- %(point_opacity)s

coefs <- coef(lm(%(yseries)s~%(xseries)s, data = plot_data))
coefs
fitcoefs = coef(lm(%(yseries)s~0 + %(xseries)s, data = plot_data))
fitlmv_Predicted <- as.numeric(fitcoefs[1])

lmv_intercept <- as.numeric(coefs[1])
lmv_Predicted <- as.numeric(coefs[2])

lm(plot_data$%(yseries)s~plot_data$%(xseries)s)
fitcoefs

xlabel <- expression(paste(plain("%(xseries)s ")*Delta*Delta*plain("G (kcal/mol)")))
ylabel <- expression(paste(plain("%(yseries)s ")*Delta*Delta*plain(G)))
rvalue <- cor(plot_data$%(yseries)s, plot_data$%(xseries)s)

minx <- min(plot_data$%(xseries)s)
maxx <- max(plot_data$%(xseries)s)
miny <- min(plot_data$%(yseries)s)
maxy <- max(plot_data$%(yseries)s)
xpos <- minx + ((maxx - minx) * 0.05)
ypos_cor <- maxy - ((maxy - miny) * 0.015)
ypos_mae <- maxy - ((maxy - miny) * 0.055)

%(plot_scale_line)s

p <- ggplot(data = plot_data, aes(x = %(xseries)s, y = %(yseries)s)) +%(plot_scale_argument)s %(extra_commands)s
    xlab("Experimental (kcal/mol)") +
    ylab("Predictions (energy units)") +
    ggtitle("%(title)s") +
    geom_point(aes(color = %(colorseries)s), alpha = I(opacity), shape = I(19)) +
    geom_abline(size = 0.25, intercept = lmv_intercept, slope = lmv_Predicted) +
    geom_abline(color="blue",size = 0.25, intercept = 0, slope = fitlmv_Predicted) +
    geom_text(hjust=0, size=4, aes(xpos, ypos_cor, fontface="plain", family = "sans", label=sprintf("R = %%0.3f", round(rvalue, digits = 4)))) +
    geom_text(hjust=0, size=4, aes(xpos, ypos_mae, fontface="plain", family = "sans", label=sprintf("MAE = %%0.3f", round(%(mae)s, digits = 4))))
p

''' % locals()



    def scatterplot_exposure(self, data, title, csv_filename):
        '''Scatterplot by exposure class.'''
        lines = ['Experimental,Predicted,Exposure,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['WildTypeExposure']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    values = c( "None" = '#777777', "B" = '%(brown)s', "E" = '%(purple)s'),
    labels = c( "None" = "N/A", "B" = "Buried", "E" = "Exposed"))''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "Exposure", title = title, plot_scale = plot_scale)


    def scatterplot_charges(self, data, title, csv_filename):
        '''Scatterplot by residue charge.'''
        lines = ['Experimental,Predicted,ResidueCharge,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['ResidueCharges']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    values = c( "None" = '#777777', "Change" = '%(cornflower_blue)s', "Polar/Charged" = 'magenta', "Hydrophobic/Non-polar" = 'green'),
    labels = c( "None" = "N/A", "Change" = "Change", "Polar/Charged" = "Polar/Charged", "Hydrophobic/Non-polar" = "Hydrophobic/Non-polar"))''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "ResidueCharge", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_volume(self, data, title, csv_filename):
        '''Scatterplot by change in volume upon mutation.'''
        lines = ['Experimental,Predicted,VolumeChange,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['VolumeChange']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    values = c( "None" = '#777777', "SL" = '%(brown)s', "LS" = '%(purple)s', 'XX' = "%(cornflower_blue)s"),
    labels = c( "None" = "N/A", "SL" = "Increase", "LS" = "Decrease", "XX" = "No change"))''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "VolumeChange", title = title, plot_scale = plot_scale)


    def scatterplot_ss(self, data, title, csv_filename):
        '''Scatterplot by secondary structure.'''
        lines = ['Experimental,Predicted,WTSecondaryStructure,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['WildTypeDSSPSimpleSSType']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    name="Secondary structure",
    values = c( "None" = '#777777', "H" = 'magenta', "S" = 'orange', "O" = '%(cornflower_blue)s'),
    labels = c( "None" = "N/A", "H" = "Helix", "S" = "Sheet", "O" = "Other"))''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "WTSecondaryStructure", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_scop_class(self, data, title, csv_filename):
        '''Scatterplot by SCOPe class.'''
        lines = ['Experimental,Predicted,WildTypeSCOPClass,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['WildTypeSCOPClass']))
        write_file(csv_filename, '\n'.join(lines))
        return self.scatterplot_color_by_series(data, colorseries = "WildTypeSCOPClass", title = title, point_opacity = 0.6)


    def scatterplot_scop_fold(self, data, title, csv_filename):
        '''Scatterplot by SCOPe fold.'''
        lines = ['Experimental,Predicted,WildTypeSCOPFold,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['WildTypeSCOPFold']))
        write_file(csv_filename, '\n'.join(lines))
        return self.scatterplot_color_by_series(data, colorseries = "WildTypeSCOPFold", title = title, point_opacity = 0.6)


    def scatterplot_scop_classification(self, data, title, csv_filename):
        '''Scatterplot by SCOPe classification.'''
        lines = ['Experimental,Predicted,WildTypeSCOPClassification,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['WildTypeSCOPClassification']))
        write_file(csv_filename, '\n'.join(lines))
        return self.scatterplot_color_by_series(data, colorseries = "WildTypeSCOPClassification", title = title, point_opacity = 0.6)


    def scatterplot_wildtype_aa(self, data, title, csv_filename):
        '''Scatterplot by wildtype residue.'''
        lines = ['Experimental,Predicted,WildTypeAA,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['WildTypeAA']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    name="Residue",
    values = c( "None" = '#808080', "A" = '#FF0000', "C" = '#BFBF00', "D" = '#008000', "E" = "#80FFFF", "F" = "#8080FF", "G" = "#BF40BF", "H" = "#A0A424", "I" = "#411BEA", "K" = "#1EAC41", "L" = "#F0C80E", "M" = "#B430E5", "N" = "#ED7651", "P" = "#19CB97", "Q" = "#362698", "R" = "#7E7EB8", "S" = "#603000", "T" = "#A71818", "V" = "#DF8020", "W" = "#E75858", "Y" = "#082008"),
    labels = c( "None" = "N/A", "A" = "A", "C" = "C", "D" = "D", "E" = "E", "F" = "F", "G" = "G", "H" = "H", "I" = "I", "K" = "K", "L" = "L", "M" = "M", "N" = "N", "P" = "P", "Q" = "Q", "R" = "R", "S" = "S", "T" = "T", "V" = "V", "W" = "W", "Y" = "Y"))
    ''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "WildTypeAA", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_mutant_aa(self, data, title, csv_filename):
        '''Scatterplot by mutant residue.'''
        lines = ['Experimental,Predicted,MutantAA,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['MutantAA']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    name="Residue",
    values = c( "None" = '#808080', "A" = '#FF0000', "C" = '#BFBF00', "D" = '#008000', "E" = "#80FFFF", "F" = "#8080FF", "G" = "#BF40BF", "H" = "#A0A424", "I" = "#411BEA", "K" = "#1EAC41", "L" = "#F0C80E", "M" = "#B430E5", "N" = "#ED7651", "P" = "#19CB97", "Q" = "#362698", "R" = "#7E7EB8", "S" = "#603000", "T" = "#A71818", "V" = "#DF8020", "W" = "#E75858", "Y" = "#082008"),
    labels = c( "None" = "N/A", "A" = "A", "C" = "C", "D" = "D", "E" = "E", "F" = "F", "G" = "G", "H" = "H", "I" = "I", "K" = "K", "L" = "L", "M" = "M", "N" = "N", "P" = "P", "Q" = "Q", "R" = "R", "S" = "S", "T" = "T", "V" = "V", "W" = "W", "Y" = "Y"))
    ''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "MutantAA", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_GP(self, data, title, csv_filename):

        lines = ['Experimental,Predicted,GP,Opacity']
        for record in data:
            residues = set([record['WildTypeAA'], record['MutantAA']])
            if None in residues:
                case_type, opacity = None, 0.5
            elif 'G' in residues or 'P' in residues:
                case_type, opacity = 'GP', 0.9
            else:
                case_type, opacity = 'Other', 0.5
            lines.append('{0},{1},{2},{3}'.format(record['Experimental'], record['Predicted'], case_type, opacity))
        write_file(csv_filename, '\n'.join(lines))

        plot_scale = '''plot_scale <- scale_color_manual(
    name="Glycine/Proline",
    values = c( "None" = '#777777', "GP" = '%(neon_green)s', "Other" = '#440077'),
    labels = c( "None" = "N/A", "GP" = "GP", "Other" = "Other"))''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "GP", title = title, plot_scale = plot_scale, point_opacity = 0.75)


    def scatterplot_pdb_res_binned(self, data, title, csv_filename):
        '''Scatterplot by binned PDB resolution.'''
        lines = ['Experimental,Predicted,PDBResolutionBin,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['PDBResolutionBin']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    name = "Resolution",
    values = c( "N/A" = '#777777', "<1.5" = '#0052aE', "1.5-2.0" = '#554C54', '2.0-2.5' = "#FFA17F", '>=2.5' = "#ce4200")
    )''' % plot_colors
        return self.scatterplot_color_by_series(data, colorseries = "PDBResolutionBin", title = title, plot_scale = plot_scale, point_opacity = 0.75)


    def scatterplot_chain_length(self, data, title, csv_filename):
        '''Scatterplot by chain length.'''
        lines = ['Experimental,Predicted,Residues,Opacity']
        for record in data:
            lines.append('{0},{1},{2},0.4'.format(record['Experimental'], record['Predicted'], record['MonomerLength']))
        write_file(csv_filename, '\n'.join(lines))
        plot_scale = '''
plot_scale <- scale_color_manual(
    name = "Resolution",
    values = c( "N/A" = '#777777', "<1.5" = '#0052aE', "1.5-2.0" = '#554C54', '2.0-2.5' = "#FFA17F", '>=2.5' = "#ce4200")
    )''' % plot_colors
        extra_commands ='\n    scale_colour_gradient(low="yellow", high="#880000") +'
        return self.scatterplot_color_by_series(data, colorseries = "Residues", title = title, plot_scale = '', point_opacity = 0.75, extra_commands = extra_commands)


if __name__ == '__main__':
    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)

    try:
        bm = BenchmarkManager(arguments)
        bm.analyze()
    except Exception, e:
        colortext.error(str(e))
        print(traceback.format_exc())

