# This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

import os
import math
import time
import tempfile
import subprocess
import traceback
import inspect
import gzip
import numpy
from scipy.stats import pearsonr, spearmanr, normaltest, ks_2samp, kstest, norm


def _id(x): pass
delete_file = [os.remove, _id][0] # for debugging, set 0 to 1
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))


def open_temp_file(path, ftype = 'w', suffix = '', prefix = ''):
    F, fname = tempfile.mkstemp(dir = path, suffix = suffix, prefix = prefix)
    output_handle = os.fdopen(F, ftype)
    return output_handle, fname


def write_temp_file(path, contents, ftype = 'w', suffix = '', prefix = ''):
    output_handle, fname = open_temp_file(path, ftype = ftype, suffix = suffix, prefix = prefix)
    output_handle.write(contents)
    output_handle.close()
    return fname


def read_file(filepath, binary = False):
    if binary:
        output_handle = open(filepath, 'rb')
    elif filepath.endswith('.gz'):
        output_handle = gzip.open(filepath, 'r')
    else:
        output_handle = open(filepath, 'r')
    contents = output_handle.read()
    output_handle.close()
    return contents


class RInterface(object):

    @staticmethod
    def _runRScript(RScript):
        rscriptname = write_temp_file(".", RScript)
        #p = subprocess.Popen(["/opt/R-2.15.1/bin/R","CMD", "BATCH", rscriptname])
        p = subprocess.Popen(["R", "CMD", "BATCH", rscriptname])
        while True:
            time.sleep(0.3)
            errcode = p.poll()
            if errcode != None:
                break
        rout = "%s.Rout" % rscriptname
        delete_file(rscriptname)

        rout_contents = None
        if os.path.exists(rout):
            rout_contents = read_file(rout)
            delete_file(rout)
        if errcode != 0:
            print(rout_contents )
            raise Exception("The R script failed with error code %d." % errcode)
        return rout_contents


    @staticmethod
    def correlation_coefficient_gplot(inputfname, output_filename, filetype, experiment_field = "Experimental"):
        '''File suffix: pearsons_r_gplot
           Description: Pearson's r
           Filename: ggplot_pearsons.R
           Priority: 1
           '''
        title = "" #Pearson's r"
        RScript = read_file(os.path.join(script_path, "ggplot_pearsons.R")) % vars()
        return RInterface._runRScript(RScript)


def create_csv(analysis_table):
    contents = '\n'.join(['ID,Experimental,Predicted'] + ['%s,%s,%s' % (str(l['ID']), str(l['Experimental']), str(l['Predicted'])) for l in analysis_table])
    return write_temp_file('.', contents)


def plot(analysis_table, output_filename, RFunction):
    #R_return_values = {}
    filetype = os.path.splitext(output_filename)[1].lower()
    assert(filetype == '.png' or filetype == '.pdf' or filetype == '.eps')
    filetype = filetype[1:]
    if len(analysis_table) <= 1:
        raise Exception("The analysis table must have at least two points.")
    else:
        input_filename = create_csv(analysis_table)
        print(input_filename )
        print(read_file(input_filename))
        try:
            R_output = RFunction(input_filename, output_filename, filetype)
            #R_return_values = RUtilities.parse_R_output(R_output)
            #for k, v in sorted(R_return_values.iteritems()):
            #    print("  %s: %s" % (str(k), str(v)))
        except Exception, e:
            print(traceback.format_exc())
            delete_file(input_filename)
            raise Exception(e)
        delete_file(input_filename)
    return output_filename


def get_ranks(values):
    '''
    Converts raw values into ranks for rank correlation coefficients
    :param values: list of values (int/float)
    :return: a dict mapping value -> rank
    '''
    ranks = {}
    sorted_values = sorted(values)
    for i in range(len(sorted_values)):
        value = sorted_values[i]
        if value not in ranks:
            ranks[value] = i + 1
    return ranks


def gamma(ranks_list1, ranks_list2):
    '''
    Goodman and Kruskal's gamma correlation coefficient
    :param ranks_list1: a list of ranks (integers)
    :param ranks_list2: a second list of ranks (integers) of equal length with corresponding entries
    :return: Gamma correlation coefficient (rank correlation ignoring ties)
    '''
    num_concordant_pairs = 0
    num_discordant_pairs = 0
    num_tied_x = 0
    num_tied_y = 0
    num_tied_xy = 0
    num_items = len(ranks_list1)
    for i in range(num_items):
        rank_1 = ranks_list1[i]
        rank_2 = ranks_list2[i]
        for j in range(i + 1, num_items):
            diff1 = ranks_list1[j] - rank_1
            diff2 = ranks_list2[j] - rank_2
            if (diff1 > 0 and diff2 > 0) or (diff1 < 0 and diff2 < 0):
                num_concordant_pairs += 1
            elif (diff1 > 0 and diff2 < 0) or (diff1 < 0 and diff2 > 0):
                num_discordant_pairs += 1
            elif diff1 == 0 and diff2 == 0:
                num_tied_xy += 1
            elif diff1 == 0:
                num_tied_x += 1
            elif diff2 == 0:
                num_tied_y += 1
    try:
        gamma_corr_coeff = float(num_concordant_pairs - num_discordant_pairs)/float(num_concordant_pairs + num_discordant_pairs)
    except:
        gamma_corr_coeff = 'n/a'
    return [num_tied_x, num_tied_y, num_tied_xy, gamma_corr_coeff]


def gamma_CC(values_list1, values_list2):
    '''
    Goodman and Kruskal's gamma correlation coefficient wrapper
    :param values_list1: a list of values
    :param values_list2: a second list of values of equal length with corresponding entries
    :return: Gamma correlation coefficient (rank correlation ignoring ties)
    '''
    ranks1 = get_ranks(values_list1)
    ranks_list1 = []
    for value in values_list1:
        rank = ranks1[value]
        ranks_list1.append(rank)
    ranks2 = get_ranks(values_list2)
    ranks_list2 = []
    for value in values_list2:
        rank = ranks2[value]
        ranks_list2.append(rank)
    return gamma(ranks_list1, ranks_list2)[3]


def fraction_correct(x_values, y_values, x_cutoff = 1.0, y_cutoff = 1.0):
    '''
    An approximation to the metric used in the Kellogg et al. paper: "The fraction correct is defined as the number of mutations categorized correctly divided by the total number of mutations in the benchmark set."
    '''
    num_points = len(x_values)
    assert(num_points == len(y_values))
    correct = 0.0
    for i in range(num_points):
        x = x_values[i]
        y = y_values[i]
        if (x >= x_cutoff) and (y >= y_cutoff): # both positive
            correct += 1.0
        elif (x <= -x_cutoff) and (y <= -y_cutoff): # both negative
            correct += 1.0
        elif (-x_cutoff < x < x_cutoff) and (-y_cutoff < y < y_cutoff): # both neutral
            correct += 1.0
    return correct / float(num_points)


def fraction_correct_fuzzy_linear_create_vector(z, z_cutoff, z_fuzzy_range):
    '''A helper function for fraction_correct_fuzzy_linear.'''
    assert(z_fuzzy_range * 2 < z_cutoff)
    if (z >= z_cutoff + z_fuzzy_range): # positive e.g. z >= 1.1
        return [0, 0, 1]
    elif (z <= -z_cutoff - z_fuzzy_range):  # negative e.g. z <= -1.1
        return [1, 0, 0]
    elif (-z_cutoff + z_fuzzy_range <= z <= z_cutoff - z_fuzzy_range): # neutral e.g. -0.9 <= z <= 0.9
        return [0, 1, 0]
    elif (-z_cutoff - z_fuzzy_range < z < -z_cutoff + z_fuzzy_range): # negative/neutral e.g. -1.1 < z < 0.9
        neutrality = (z + z_cutoff + z_fuzzy_range) / (z_fuzzy_range * 2)
        zvec = [1 - neutrality, neutrality, 0]
    elif (z_cutoff - z_fuzzy_range < z < z_cutoff + z_fuzzy_range): # neutral/positive e.g. 0.9 < z < 1.1
        positivity = (z - z_cutoff + z_fuzzy_range) / (z_fuzzy_range * 2)
        zvec = [0, 1 - positivity, positivity]
    else:
        raise Exception('Logical error.')
    # normalize the vector
    length = math.sqrt(numpy.dot(zvec, zvec))
    return numpy.divide(zvec, length)


def fraction_correct_fuzzy_linear(x_values, y_values, x_cutoff = 1.0, x_fuzzy_range = 0.1, y_scalar = 1.0):
    '''
    A version of fraction_correct which is more forgiving at the boundary positions.
    In fraction_correct, if the x value is 1.01 and the y value is 0.99 (with cutoffs of 1.0) then that pair evaluates
    to zero despite the results being very close to each other.
    This function corrects for the boundary by overlapping the ranges and attenuating the endpoints.
    This version of the function uses a linear approach - a classification (positive, negative, neutral resp. P, N, X)
    is 1 for some range of values, 0 for a separate range, and between 0 and 1 for the in-between range i.e.
        N       X        P
       ----\  /----\  /-----
            \/      \/
            /\      /\
       ----/  \----/  \----
    '''
    num_points = len(x_values)
    assert(num_points == len(y_values))
    correct = 0.0
    for i in range(num_points):
        x = x_values[i]
        y = y_values[i]
        y_cutoff = x_cutoff * y_scalar
        y_fuzzy_range = x_fuzzy_range * y_scalar
        xvec = fraction_correct_fuzzy_linear_create_vector(x, x_cutoff, x_fuzzy_range)
        yvec = fraction_correct_fuzzy_linear_create_vector(y, y_cutoff, y_fuzzy_range)
        correct += numpy.dot(xvec, yvec)
    return correct / float(num_points)


def mae(x_values, y_values):
    '''Mean absolute/unsigned error.'''
    num_points = len(x_values)
    assert(num_points == len(y_values) and num_points > 0)
    return numpy.sum(numpy.apply_along_axis(numpy.abs, 0, numpy.subtract(x_values, y_values))) / float(num_points)


def get_xy_dataset_statistics(x_values, y_values, fcorrect_x_cutoff = 1.0, fcorrect_y_cutoff = 1.0, x_fuzzy_range = 0.1, y_scalar = 1.0):
    '''A function which takes two lists of values of equal length with corresponding entries and returns a dict containing
       a variety of metrics.'''
    assert(len(x_values) == len(y_values))
    return dict(
        pearsonr = pearsonr(x_values, y_values),
        spearmanr = spearmanr(x_values, y_values),
        gamma_CC = gamma_CC(x_values, y_values),
        MAE = mae(x_values, y_values),
        normaltestx = normaltest(x_values),
        normaltesty = normaltest(y_values),
        kstestx = kstest(x_values, 'norm'),
        kstesty = kstest(y_values, 'norm'),
        ks_2samp = ks_2samp(x_values, y_values),
        fraction_correct = fraction_correct(x_values, y_values, x_cutoff = fcorrect_x_cutoff, y_cutoff = fcorrect_y_cutoff),
        fraction_correct_fuzzy_linear = fraction_correct_fuzzy_linear(x_values, y_values, x_cutoff = fcorrect_x_cutoff, x_fuzzy_range = x_fuzzy_range, y_scalar = y_scalar),
    )


if __name__ == '__main__':
    import random

    # Create a random dataset that is skewed towards correlation
    test_table = []
    for x in range(1, 1001):
        exp = random.uniform(-5, 5)
        if random.random() < 0.2:
            pred = exp + random.uniform(-1, 1)
            test_table.append(dict(ID = x, Experimental = exp, Predicted = pred))
        elif random.random() < 0.5:
            pred = exp + random.uniform(-2, 2)
            test_table.append(dict(ID = x, Experimental = exp, Predicted = pred))
        else:
            test_table.append(dict(ID = x, Experimental = exp, Predicted = random.uniform(-5, 5)))

    plot(test_table, 'test.png', RInterface.correlation_coefficient_gplot)