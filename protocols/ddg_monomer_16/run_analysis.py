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
    run_analysis.py [options]...

Options:

    -o --output_directory OUTPUT_DIR
        The path to a directory previously generated from the run_preminimization script. This defaults to the most recent directory in job_output, if this exists.

    -p --scatterplot_filename SCATTERPLOT_FILE
        The filename of the scatterplot to be generated in the output directory (unless --skip_analysis is set). [default: scatterplot.png] 
    
    --force
        When this option is set, the most recent directory in job_output, if it exists, will be used without prompting the user.

    --skip_analysis
        When this option is set, the analysis script is not invoked once the analysis files are generated.

    --take_lowest N
        When this option is set, the average of the N lowest-scoring (most stable) mutant and wildtypes structures are used to calculate the DDG value. [default: 3]

    --use_single_reported_value
        By default, the analysis takes the three lowest-scoring mutant structures and the three lowest-scoring wildtype structures to calculate the DDG value. This approach was taken by Kellogg et al. and should reduce stochastic noise. If this option is set, the single value reported by ddg_monomer is used instead. We do not recommend using this option.

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
import glob
import cPickle as pickle
import getpass
import gzip
from rosetta.write_run_file import process as write_run_file
from analysis.libraries import docopt
from analysis.stats import read_file, read_file_lines, write_file, prompt_yn
from run_ddg import task_subfolder as ddg_task_subfolder
try:
    import json
except:
    import simplejson as json


generated_csv = 'benchmark_output.csv'
generated_json = 'benchmark_output.json'
ddGMover_regex = re.compile("^protocols.moves.ddGMover:\s*mutate\s*.*?\s*wildtype_dG\s*is:\s*.*?and\s*mutant_dG\s*is:\s*.*?\s*ddG\s*is:\s*(.*)$")

van_der_Waals_volumes = dict(
    G = 48,
    A = 67,
    S = 73,
    C = 86,
    P = 90,
    D = 91,
    T = 93,
    N = 96,
    V = 105,
    E = 109,
    Q = 114,
    H = 118,
    I = 124,
    L = 124,
    M = 124,
    F = 135,
    K = 135,
    Y = 141,
    R = 148,
    W = 163,
)

def read_stdout(ddg_output_path):
    try:
        output_file = os.path.join(ddg_output_path, 'rosetta.out.gz')
        rosetta_output = read_file_lines(output_file)
    except:
        raise Exception('An error occurred reading the output file %s.' % output_file)
    return rosetta_output


def determine_SL_class(record):
    '''Returns:
        - SL if all of the mutations are from smaller (by volume) residues to larger residues
        - LS if all of the mutations are from larger residues to smaller residues
        - XX if all of the mutations are to and from residues of the same size e.g. M -> L
        - O  if none of the cases above apply (this occurs when records with multiple mutations mix the types above).'''

    mutation_types = set()
    for m in record['Mutations']:
        wt_vol, mut_vol = van_der_Waals_volumes[m['WildTypeAA']], van_der_Waals_volumes[m['MutantAA']]
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


def extract_predicted_ddg(ddg_output_path):
    rosetta_output_lines = read_stdout(ddg_output_path)
    ddGMover_lines = [l for l in rosetta_output_lines if l.strip() and l.startswith("protocols.moves.ddGMover: mutate")]
    assert(len(ddGMover_lines) == 1)
    ddGMover_line = ddGMover_lines[0].strip()
    mtchs = ddGMover_regex.match(ddGMover_line)
    assert(mtchs)
    return float(mtchs.group(1))


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


def get_ddg_monomer_scores_per_structure(ddg_output_path):
    '''Returns a dict mapping the DDG scores from a ddg_monomer run to a list of structure numbers.'''

    rosetta_output_lines = read_stdout(ddg_output_path)

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


def extract_data(ddg_output_path, take_lowest):
    '''Note: This function is written to handle data from Rosetta at the time of writing. This will need to be altered to
       work with certain older revisions and may need to be adapted to work with future revisions.'''
    scores = {}
    scores['DDG'] = extract_predicted_ddg(ddg_output_path)
    scores['DDG_components'] = extract_summary_data(ddg_output_path)
    for k, v in get_ddg_monomer_scores_per_structure(ddg_output_path).iteritems():
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


if __name__ == '__main__':
    import pprint
    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)

    take_lowest = None
    try:
        assert(arguments['--take_lowest'][0].isdigit())
        take_lowest = int(arguments['--take_lowest'][0])
        assert(take_lowest)
    except:
        print('take_lowest must be a positive non-zero integer: "%s" was passed.' % arguments['--take_lowest'][0])
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
                    print('\nRunning analysis in %s.' % most_recent_directory)
                else:
                    answer = prompt_yn('\nNo output path was specified. Use %s (y/n)?' % most_recent_directory)
                if not answer:
                    print('No output path was specified. Exiting.\n')
                    sys.exit(1)
                output_dir = most_recent_directory
            else:
                print('No output could be found in the job_output directory. Exiting.\n')
                sys.exit(1)

    # Set the job output directories
    output_data_dir = os.path.join(output_dir, 'data')
    ddg_data_dir = os.path.join(output_dir, ddg_task_subfolder)
    for p in [output_data_dir, ddg_data_dir]:
        if not os.path.exists(p):
            raise Exception('The folder %s was expected to exist after the preminimization step but could not be found.' % p)

    # Read in the dataset file
    try:
        dataset_filepath = os.path.join(output_dir, 'dataset.json')
        dataset = json.loads(read_file(dataset_filepath))
        dataset_cases = {} 
        for dsc in dataset['data']:
            dataset_cases[dsc['RecordID']] = dsc
    except Exception, e:
        raise Exception('An error occurred parsing the JSON file: %s..' % str(e))

    # Save the benchmark data to file
    benchmark_data_filepath = os.path.join(output_dir, 'benchmark_data.json')
    print('Creating %s which contains component and summary scores for each case and generated structure.' % benchmark_data_filepath)
    analysis_data = {}
    job_dirs = glob.glob(os.path.join(ddg_data_dir, '*'))
    for jd in job_dirs:
        if os.path.isdir(jd):
            jdirname = os.path.split(jd)[1]
            if jdirname.isdigit():
                record_id = int(jdirname)
                analysis_data[record_id] = extract_data(jd, take_lowest)
    write_file(benchmark_data_filepath, json.dumps(analysis_data, indent = 4))

    # Create XY data
    analysis_csv_input_filepath = os.path.join(output_dir, 'analysis_input.csv')
    analysis_json_input_filepath = os.path.join(output_dir, 'analysis_input.json')
    print('Creating input files %s and %s for the analysis script.' % (analysis_csv_input_filepath, analysis_json_input_filepath))
    if len(analysis_data) > len(dataset_cases):
        raise Exception('ERROR: There seems to be an error - there are more predictions than cases in the dataset. Exiting.')
    elif len(analysis_data) < len(dataset_cases):
        print('\nWARNING: %d cases missing for analysis; there are %d predictions in the output directory but %d cases in the dataset. The analysis below does not cover the complete dataset.\n' % (len(dataset_cases) - len(analysis_data), len(analysis_data), len(dataset_cases)))

    # ddg_analysis_type can be set to 'DDG' or 'DDG_Top3'.
    # 'DDG' uses the value reported by ddg_monomer.
    # 'DDG_Top3' (generated by default) uses the metric from Kellogg et al. based on the three lowest scoring mutant structures and the three lowest scoring wildtype structures
    ddg_analysis_type = 'DDG_Top%d' % take_lowest
    if arguments['--use_single_reported_value']:
        ddg_analysis_type = 'DDG'
        print('\nThe predicted DDG value per case is the single DDG value reported by the ddg_monomer application.')
    elif take_lowest == 3:
        print('\nThe predicted DDG value per case is computed using the %d lowest-scoring mutant structures and the %d lowest-scoring wildtype structures as in the paper by Kellogg et al.' % (take_lowest, take_lowest))
    else:
        print('\nThe predicted DDG value per case is computed using the %d lowest-scoring mutant structures and the %d lowest-scoring wildtype structures.' % (take_lowest, take_lowest))


    json_records = dict(
        All = [],
            ByVolume = dict(
                SL = [], LS = [], XX = [], O = []
            ),
            GP = {
                True : [], False : []
            },
    )
    csv_file = ['#Experimental,Predicted,ID']

    for record_id, predicted_data in sorted(analysis_data.iteritems()):
        record = dataset_cases[record_id]
        json_record = dict(Experimental = record['DDG'], Predicted = predicted_data[ddg_analysis_type], ID = record_id)
        json_records['ByVolume'][determine_SL_class(record)].append(json_record)
        json_records['GP'][has_G_or_P(record)].append(json_record)
        json_records['All'].append(json_record)
        csv_file.append('%s,%s,%s' % (str(record['DDG']), str(predicted_data[ddg_analysis_type]), str(record_id)))
    write_file(analysis_csv_input_filepath, '\n'.join(csv_file))
    write_file(analysis_json_input_filepath, json.dumps(json_records['All'], indent = 4))

    by_volume_descriptions = dict(
        SL = 'small-to-large mutations',
        LS = 'large-to-small mutations',
        XX = 'no change in volume',
    )

    if not arguments['--skip_analysis']:
        from analysis.stats import get_xy_dataset_statistics, plot, format_stats_for_printing, RInterface
        correlation_coefficient_scatterplotplot = RInterface.correlation_coefficient_gplot

        # Set up the output filename
        print('\nRunning analysis:')


        volume_groups = {}
        for k, v in van_der_Waals_volumes.iteritems():
            volume_groups[v] = volume_groups.get(v, [])
            volume_groups[v].append(k)

        print('\n\nSection 1. Breakdown by volume.')
        print('A case is considered a small-to-large (resp. large-to-small) mutation if all of the wildtype residues have a smaller (resp. larger) van der Waals volume than the corresponding mutant residue. The order is defined as %s so some cases are considered to have no change in volume e.g. MET -> LEU.' % (' < '.join([''.join(sorted(v)) for k, v in sorted(volume_groups.iteritems())])))
        for subcase in ('XX', 'SL', 'LS'):
            print('\n' + '*'*10 + (' Statistics - %s (%d cases)' % (by_volume_descriptions[subcase], len(json_records['ByVolume'][subcase]))) +'*'*10)
            print(format_stats_for_printing(get_xy_dataset_statistics(json_records['ByVolume'][subcase])))

        print('\n\nSection 2. Separating out mutations involving glycine or proline.')
        print('This cases may involve changes to secondary structure so we separate them out here.')
        print('\n' + '*'*10 + (' Statistics - cases with G or P (%d cases)' % (len(json_records['GP'][True]))) +'*'*10)
        print(format_stats_for_printing(get_xy_dataset_statistics(json_records['GP'][True])))
        print('\n' + '*'*10 + (' Statistics - cases without G or P (%d cases)' % (len(json_records['GP'][False]))) +'*'*10)
        print(format_stats_for_printing(get_xy_dataset_statistics(json_records['GP'][False])))

        print('\n\nSection 3. Entire dataset.')
        print('\n' + '*'*10 + (' Statistics - complete dataset (%d cases)' % len(json_records['All'])) +'*'*10)
        print(format_stats_for_printing(get_xy_dataset_statistics(json_records['All'])))

        output_filename = os.path.join(output_dir, arguments['--scatterplot_filename'][0])
        print('\nSaving scatterplot to %s.\n' % output_filename)
        plot(json_records['All'], output_filename, correlation_coefficient_scatterplotplot)
        

    print('')
