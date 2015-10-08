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
A simple application to outputs statistics and a scatterplot for the given XY data set.
See reports/ddg_analyzer.py for a more fully-featured application for detailed analyses and
protocols/ddg_monomer_16/run_analysis.py for an example of how to specialize this application to a particular computational
method.

Usage:
    analyze.py [options] <inputfile>...

Options:
    -o FILE --output FILE
        File name for the generated scatterplot [default: scatterplot.png]

The input file should be a comma-separated values file where the first three columns are:
  ID,Experimental,Predicted

Authors:
    Shane O'Connor
"""

import sys
import os

from libraries import docopt

from klab.stats import get_xy_dataset_statistics, format_stats_for_printing
from klab.fs.fsio import read_file
from klab.benchmarking.analysis.plot import plot
from klab.plot.rtools import RInterface


correlation_coefficient_scatterplotplot = RInterface.correlation_coefficient_gplot


def read_json(filename):
    try:
        try:
            import json
        except:
            import simplejson as json
        return json.loads(read_file(filename))
    except:
        return None


def parse_csv(filename):
    separator = ','
    if filename.endswith('.tsv'):
        separator = '\t'
    try:
        table = []
        id = 1
        contents = read_file(filename)
        lines = [l.strip().split(separator) for l in contents.split('\n') if l.strip() and not(l.strip().startswith('#'))]
        for linetokens in lines:
            if len(linetokens) >= 3:
                table.append(dict(Experimental = float(linetokens[0]), Predicted = float(linetokens[1]), ID = str(linetokens[2])))
            elif len(linetokens) == 2:
                table.append(dict(Experimental = float(linetokens[0]), Predicted = float(linetokens[1]), ID = id))
                id += 1
            else:
                raise Exception('At least two columns (experimental DDG, predicted DDG) are expected.')
        return table
    except Exception, e:
        raise Exception('An exception occurred parsing the CSV/TSV file: %s' % str(e))


if __name__ == '__main__':
    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)

    # Read file input file
    input_filename = arguments['<inputfile>'][0]
    if not os.path.exists(input_filename):
        print('Error: the input file %s does not exist.' % input_filename)
        sys.exit(2)
    analysis_table = read_json(input_filename)
    if not analysis_table:
        analysis_table = parse_csv(input_filename)

    # Set up the output filename
    output_filename = arguments['--output']
    output_filename_ext = os.path.splitext(output_filename)[1].lower()
    if output_filename_ext not in ['.png', '.pdf']: # todo: check eps output ('.eps')
        output_filename += '.png'

    print('\n' + '*'*10 + ' Statistics ' +'*'*10)
    print(format_stats_for_printing(get_xy_dataset_statistics(analysis_table)))

    print('\nSaving scatterplot to %s.\n' % output_filename)
    plot(analysis_table, output_filename, correlation_coefficient_scatterplotplot)

