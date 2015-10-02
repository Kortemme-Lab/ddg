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

"""\
The script compiles the analysis files for the benchmark run and then invokes the analysis script.

 Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein
 structure and stability. 2011. Proteins. 79(3):830-8. doi: 10.1002/prot.22921.

Usage:
    run_analysis.py <analysis_directory> [-b BENCHMARK_RUN_DIRECTORY...] [-n BENCHMARK_RUN_NAME...] [options]

Options:

    -j --benchmark_run_json BENCHMARK_RUN_JSON
        The path to a JSON file describing the benchmark runs being analyzed. Either this option or the -b and -n options
        should be used but not both. If this option is used then more information will be added to the title slide of the
        generated report. An example of the file format is provided in benchmark_runs.json.example.

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
        When this option is set, the benchmark_data.json and the pandas HDF5 file is not regenerated. This saves time on
        subsequent calls to this analysis script but is disabled by default.

    --force
        When this option is set, the most recent directory in job_output, if it exists, will be used without prompting the user.

    --recreate_graphs
        When this option is set, all plots are regenerated. This is disabled by default to allow subsequent analysis to
        run faster.

    -G --do_not_generate_plots
        When this option is set, the graphical plots are not generated.

    -R --do_not_report_analysis
        When this option is set, the analyses are not printed to screen.

    --silent
        When this option is set, nothing is written to the terminal.

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
"""

import sys
import os
import re
import glob
import traceback
import pandas
try: import json
except: import simplejson as json

from analysis.libraries import docopt
from analysis.libraries import colortext
from analysis.stats import read_file, read_file_lines, write_file, prompt_yn
from analysis.libraries.loggers import ReportingObject
from run_ddg import task_subfolder as ddg_task_subfolder


# This module contains two main classes.
#
# The DDGBenchmarkManager class reads in the script options, extracts the score data from the benchmark runs, and controls
# the analysis. The score extraction is not implemented by this class.
#
# DDGMonomerBenchmarkManager is a concrete implementation of DDGBenchmarkManager with the score extraction code for the
# ddg_monomer Rosetta application.


class DDGBenchmarkManager(ReportingObject):
    '''This class is responsible for extract the data from benchmark runs and creating BenchmarkRun objects which can then be called to perform analysis.'''


    def __init__(self, arguments):
        '''Read the command-line options and extract the raw data from the benchmark runs.'''

        # Declare object variables
        self.benchmark_run_data = {}            # a mapping from benchmark names to BenchmarkRun objects
        self.analysis_directory = None
        self.take_lowest = None
        self.burial_cutoff = None
        self.include_derived_mutations = None
        self.stability_classication_x_cutoff, self.stability_classication_y_cutoff = None, None

        # Parse command-line arguments and extract/load the benchmark input data
        self.parse_arguments(arguments)
        self.setup()


    def parse_arguments(self, arguments):
        '''Read command-line options.'''

        # Output directory. The path where the analysis files should be created.
        self.analysis_directory = os.path.abspath(arguments['<analysis_directory>'])

        # Create analysis directory
        if not(os.path.exists(self.analysis_directory)):
            try:
                os.mkdir(self.analysis_directory)
                assert(os.path.exists(self.analysis_directory))
            except Exception, e:
                raise colortext.Exception('An exception occurred creating the directory %s.' % self.analysis_directory)

        # Plot-generation option
        self.generate_plots = not(arguments['--do_not_generate_plots'])
        self.report_analysis = not(arguments['--do_not_report_analysis'])
        self.silent = arguments['--silent']

        # Whether or not we include records marked as derived in the analysis
        self.include_derived_mutations = arguments['--include_derived_mutations']

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
        benchmark_run_json = arguments['--benchmark_run_json']
        benchmark_run_descriptions = []
        benchmark_run_credits = []

        if benchmark_run_json and (len(benchmark_run_directories) > 0 and len(benchmark_run_directories) > 0):
            raise colortext.Exception('The -j (--benchmark_run_json) option cannot be used if either of the -b (--benchmark_run_directory) or -n (--benchmark_run_name) options are used.')

        if benchmark_run_json:
            benchmark_json_runs = []
            if not os.path.exists(benchmark_run_json):
                raise colortext.Exception('The file {0} does not exist.'.format(benchmark_run_json))
            try:
                benchmark_run_json = json.loads(read_file(benchmark_run_json))
            except Exception, e:
                raise colortext.Exception('An error occurred parsing {0}: "{1}".'.format(benchmark_run_json, str(e)))

            try:
                c = [benchmark_run.get('Order') for benchmark_run in benchmark_run_json]
                if None in c:
                    c.remove(None)
                c = [int(x) for x in c if str(x).isdigit()]
                if c:
                    c = max(c) + 1
                else:
                    c = 1
                print(c)
                for benchmark_run in benchmark_run_json:
                    if benchmark_run.get('Order') and str(benchmark_run.get('Order')).isdigit():
                        order = int(benchmark_run['Order'])
                    else:
                        order = c
                        c += 1
                    benchmark_json_runs.append((order, benchmark_run['Path'], benchmark_run['Name'], benchmark_run.get('Description', ''), benchmark_run.get('Credit', '')))
                benchmark_json_runs = sorted(benchmark_json_runs)
                benchmark_run_directories = [os.path.normpath(os.path.expanduser(x[1])) for x in benchmark_json_runs]
                benchmark_run_names = [x[2] for x in benchmark_json_runs]
                benchmark_run_descriptions = [x[3] for x in benchmark_json_runs]
                benchmark_run_credits = [x[4] for x in benchmark_json_runs]
            except Exception, e:
                raise colortext.Exception('An error occurred parsing {0}: "{1}".'.format(benchmark_run_json, str(e)))
        else:
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
                            self.log('\nRunning analysis from the run in %s.' % most_recent_directory)
                        else:
                            answer = prompt_yn('\nNo benchmark run path was specified. Use %s (y/n)?' % most_recent_directory)
                        if not answer:
                            raise colortext.Exception('No benchmark run was specified. Exiting.\n')
                        benchmark_run_directories = [most_recent_directory]
                        if len(benchmark_run_names) == 0:
                            benchmark_run_names.append('')
                            colortext.warning('No benchmark_run_name argument was specified. Defaulting to a blank name.') # we write to stdout here despite the user option
                        elif len(benchmark_run_names) > 1:
                            raise colortext.Exception('No benchmark_run_directory was specified (one recent directory was chosen automatically) but multiple --benchmark_run_name arguments were specified.')

        # Check for uniqueness
        if benchmark_run_directories and (len(benchmark_run_directories) != len(set(benchmark_run_directories))):
            raise colortext.Exception('The benchmark run directories must be unique.')
        if benchmark_run_names and (len(benchmark_run_names) != len(set(benchmark_run_names))):
            raise colortext.Exception('The benchmark names e options must be unique.')

        # Make sure the arrays are of the same size
        assert(len(benchmark_run_directories) == len(benchmark_run_names))
        if not benchmark_run_descriptions:
            benchmark_run_descriptions = ['' for x in range(len(benchmark_run_directories))]
        if not benchmark_run_credits:
            benchmark_run_credits = ['' for x in range(len(benchmark_run_directories))]
        assert(len(benchmark_run_directories) == len(benchmark_run_descriptions))
        assert(len(benchmark_run_directories) == len(benchmark_run_credits))

        # At least one benchmark run must be chosen (including the published benchmark run)
        if not use_published_data and len(benchmark_run_directories) == 0:
            raise colortext.Exception('No benchmark runs were specified and no output could be found in the job_output directory. Exiting.\n')

        if self.analysis_directory in benchmark_run_directories:
            raise colortext.Exception('The analysis output directory {0} cannot be one of the input benchmark directories.'.format(self.analysis_directory))

        # Load in the data for each benchmark
        for x in range(len(benchmark_run_directories)):

            benchmark_run_directory = benchmark_run_directories[x]
            benchmark_run_name = benchmark_run_names[x]
            benchmark_run_description = benchmark_run_descriptions[x]
            benchmark_run_credit = benchmark_run_credits[x]
            self.log('Setting up the analysis data for {0} in {1}.'.format(benchmark_run_name, benchmark_run_directory), fn = colortext.message)

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
                dataset_description = dataset.get('description', 'Unknown dataset')
                dataset_cases = {}
                for dsc in dataset['data']:
                    dataset_cases[dsc['RecordID']] = dsc
            except Exception, e:
                raise colortext.Exception('An error occurred parsing the JSON file: %s..' % str(e))

            # Read the previously extracted benchmark data and structural scores (benchmark_data.json) from file or create that file if it does not exist
            benchmark_data_filepath = os.path.join(benchmark_run_directory, 'benchmark_data.json')
            analysis_data = {}
            if not(os.path.exists(benchmark_data_filepath)) or arguments.get('--use_existing_benchmark_data') == 0:
                self.log('Creating %s which contains component and summary scores for each case and generated structure.' % benchmark_data_filepath, fn = colortext.warning)

                job_dirs = sorted([jd for jd in glob.glob(os.path.join(ddg_data_dir, '*')) if os.path.isdir(jd) and os.path.split(jd)[1].isdigit()])
                c, num_dirs = 0.0, float(len(job_dirs)) / 100.0
                for jd in job_dirs:
                    jdirname = os.path.split(jd)[1]
                    record_id = int(jdirname)
                    c += 1.0
                    self.log('\rProgress: {0:d}%'.format(int(c/num_dirs)), colortext.wcyan)
                    if not self.silent: sys.stdout.flush()
                    analysis_data[record_id] = self.__class__.extract_data(jd, self.take_lowest)
                self.log('\rWriting to file...', colortext.wlightpurple)
                if not self.silent: sys.stdout.flush()
                write_file(benchmark_data_filepath, json.dumps(analysis_data, indent = 4, sort_keys=True))
                self.log('\r', colortext.wgrey)
                if not self.silent: sys.stdout.flush()
            else:
                self.log('Found an existing benchmark_data.json file containing component and summary scores for each case and generated structure:', colortext.warning)
                self.log('\r...loading', colortext.wyellow)
                if not self.silent: sys.stdout.flush()
                analysis_data_ = json.loads(read_file(benchmark_data_filepath))
                for k, v in analysis_data_.iteritems():
                    analysis_data[int(k)] = v
                self.log('\r', colortext.wgrey)
                if not self.silent: sys.stdout.flush()

            self.benchmark_run_data[benchmark_run_name] = BenchmarkRun(
                benchmark_run_name,
                benchmark_run_directory,
                self.analysis_directory,
                dataset_cases,
                analysis_data,
                arguments['--use_single_reported_value'],
                description = benchmark_run_description,
                dataset_description = dataset_description,
                credit = benchmark_run_credit,
                include_derived_mutations = self.include_derived_mutations,
                take_lowest = self.take_lowest,
                generate_plots = self.generate_plots,
                report_analysis = self.report_analysis,
                silent = self.silent,
                burial_cutoff = self.burial_cutoff,
                stability_classication_x_cutoff = self.stability_classication_x_cutoff,
                stability_classication_y_cutoff = self.stability_classication_y_cutoff,
                use_existing_benchmark_data = arguments['--use_existing_benchmark_data'],
                recreate_graphs = arguments['--recreate_graphs']
            )


    def analyze(self):
        '''Runs the analysis for the different benchmark runs.'''

        benchmark_runs = sorted(self.benchmark_run_data.keys())
        # todo: add published data here

        # Create the dataframes for each benchmark
        for benchmark_run_name, br in sorted(self.benchmark_run_data.iteritems()):
            self.log('\nExtracting the run-data the analysis for {0}.'.format(benchmark_run_name), colortext.message)
            br.create_dataframe()

        # Run the individual analysis
        for benchmark_run_name, br in sorted(self.benchmark_run_data.iteritems()):
            self.log('\nRunning the analysis for {0}.'.format(benchmark_run_name), colortext.message)
            br.analyze()

        # Compare the benchmark runs against each other
        for x in range(len(benchmark_runs)):
            for y in range(x + 1, len(benchmark_runs)):
                pass
                # todo: run comparative analysis


    @staticmethod
    def extract_data(ddg_output_path, take_lowest):
        raise Exception('Abstract method.')



class DDGMonomerBenchmarkManager(DDGBenchmarkManager):


    def __init__(self):
        super(DDGMonomerBenchmarkManager, self).__init__(self, arguments)

        # Regex used to extract data from Rosetta output
        self.__class__.ddGMover_regex = re.compile("^protocols.moves.ddGMover:\s*mutate\s*.*?\s*wildtype_dG\s*is:\s*.*?and\s*mutant_dG\s*is:\s*.*?\s*ddG\s*is:\s*(.*)$")


    ###
    #  Data-extraction methods
    ###


    @staticmethod
    def extract_data(ddg_output_path, take_lowest):
        '''Extract the data for one dataset case.
           Note: This function is written to handle data from Rosetta at the time of writing. This will need to be altered to
           work with certain older revisions and may need to be adapted to work with future revisions.'''
        scores = {}
        rosetta_output_lines = DDGMonomerBenchmarkManager.read_stdout(ddg_output_path)
        scores['Errors'] = DDGMonomerBenchmarkManager.extract_errors(rosetta_output_lines)
        scores['DDG'] = DDGMonomerBenchmarkManager.extract_predicted_ddg(rosetta_output_lines)
        scores['DDG_components'] = DDGMonomerBenchmarkManager.extract_summary_data(ddg_output_path)
        for k, v in DDGMonomerBenchmarkManager.get_ddg_monomer_scores_per_structure(rosetta_output_lines).iteritems():
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
        mtchs = DDGMonomerBenchmarkManager.ddGMover_regex.match(ddGMover_line)
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




if __name__ == '__main__':

    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)

    try:
        bm = DDGMonomerBenchmarkManager(arguments)
        bm.analyze()
    except Exception, e:
        colortext.error(str(e))
        print(traceback.format_exc())

