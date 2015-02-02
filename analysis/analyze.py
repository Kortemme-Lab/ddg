#!/usr/bin/env python2
# This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

"""\
Outputs statistics and a scatterplot for the given XY data set.

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
from stats import get_xy_dataset_statistics, plot, read_file, RInterface
correlation_coefficient_scatterplotplot = RInterface.correlation_coefficient_gplot

if __name__ == '__main__':
    try:
        arguments = docopt.docopt(__doc__.format(**locals()))
    except Exception, e:
        print('Failed while parsing arguments: %s.' % str(e))
        sys.exit(1)
   
    # Read file input file
    input_filename = arguments['<inputfile>'][0]
    if not os.path.exists(input_filename):
        print('Error') # todo:
        sys.exit(2)
    analysis_table = read_file(input_filename)

    # Set up the output filename
    output_filename = arguments['--output']
    output_filename_ext = os.path.splitext(output_filename)[1].lower()
    if output_filename_ext not in ['png', 'pdf', 'eps']:
        output_filename += '.png'

    plot(analysis_table, output_filename, correlation_coefficient_scatterplotplot)

