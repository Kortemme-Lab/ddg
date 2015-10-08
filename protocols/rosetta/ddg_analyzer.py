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


import sys
import os
import glob
try: import json
except: import simplejson as json

import pandas

from tools import colortext
from tools.tui.utils import prompt_yn
from tools.fs.fsio import read_file, write_file
from tools.loggers.simple import ReportingObject
from tools.benchmarking.analysis.ddg_monomeric_stability_analysis import BenchmarkRun as DefaultBenchmarkRunClass


# This module contains one main class.
#
# The DDGBenchmarkManager class reads in the script options, extracts the score data from the benchmark runs, and controls
# the analysis.
# This is meant to be an abstract class which does not implement the score extraction code and which can be used by different
# computational methods. It depends on some conditions which the setup script must satisfy - see the comments in the setup
# function. An example of how to use this class can be found in protocols/ddg_monomer_16/run_analysis.py.


class DDGBenchmarkManager(ReportingObject):
    '''This class is responsible for extract the data from benchmark runs and creating BenchmarkRun objects which can then be called to perform analysis.'''


    def __init__(self, arguments, BenchmarkRunClass = DefaultBenchmarkRunClass):
        '''Read the command-line options and extract the raw data from the benchmark runs.'''

        # Declare object variables
        self.benchmark_run_data = {}            # a mapping from benchmark names to BenchmarkRun objects
        self.analysis_directory = None
        self.BenchmarkRunClass = BenchmarkRunClass
        self.take_lowest = None
        self.burial_cutoff = None
        self.include_derived_mutations = None
        self.stability_classication_x_cutoff, self.stability_classication_y_cutoff = None, None

        # Parse command-line arguments and extract/load the benchmark input data
        self.arguments = arguments
        self.parse_arguments()
        self.setup()


    def parse_arguments(self):
        '''Read command-line options.'''

        arguments = self.arguments

        # Output directory. The path where the analysis files should be created.
        self.analysis_directory = os.path.abspath(arguments['<analysis_directory>'])

        # Output directory. The path where the analysis files should be created.
        self.ddg_task_subfolder = arguments['--data_subfolder']

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
             - the dataset used to run the benchmark is contained in a file called dataset.json in the output directory (this should be created when the benchmark run is set up)
             - the ddg results are stored in a subfolder <ddg_task_subfolder> (a default should be by the docopt script - see protocols/ddg_monomer_16/run_analysis.py for an example)
           Postconditions set up by this script:
             - the extracted benchmark data (analysis_data) is:
                - stored in a file called benchmark_data.json in the benchmark run root directory; or
                - created in that location by this script.
           The loaded data is stored inside BenchmarkRun objects for analysis and cross-analysis. For convenience, the
           analyze method of this class performs all possible single and pair-wise analysis.
        '''

        # Determine the benchmark run directories and user-specified names
        arguments = self.arguments
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
            raise colortext.Exception('The benchmark names options must be unique.')

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
            ddg_data_dir = os.path.join(benchmark_run_directory, self.ddg_task_subfolder)
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

            self.benchmark_run_data[benchmark_run_name] = self.BenchmarkRunClass(
                benchmark_run_name,
                benchmark_run_directory,
                self.analysis_directory,
                dataset_cases,
                analysis_data,
                use_single_reported_value = arguments['--use_single_reported_value'],
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

