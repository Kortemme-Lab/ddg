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

    --data_subfolder
        The subfolder of the benchmark run directories where data is stored if this is not found in the root folder. [default: ddg]

Authors:
    Shane O'Connor
"""

import sys
import os
import re
import traceback
try: import json
except: import simplejson as json

from libraries import docopt

from klab import colortext
from klab.fs.fsio import read_file_lines

from rosetta.ddg_analyzer import DDGBenchmarkManager


# This module contains one main class.
#
# DDGMonomerBenchmarkManager is a concrete implementation of DDGBenchmarkManager with the score extraction code for the
# ddg_monomer Rosetta application.


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

