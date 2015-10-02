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


# todo: This should be moved into another directory.


import os
import shutil
import numpy
import pprint
import subprocess
import shlex
import pandas
import StringIO
try: import json
except: import simplejson as json

from analysis.libraries import colortext
from analysis.libraries.colors import rgb_colors as plot_colors
from analysis.stats import read_file, write_file, fraction_correct, fraction_correct_pandas, add_fraction_correct_values_to_dataframe, get_xy_dataset_statistics_pandas, format_stats_for_printing, RInterface, plot_pandas
from analysis.libraries.loggers import ReportingObject


# This module contains one main class.
#
# The BenchmarkRun class creates a dataframe containing the raw data used for analysis. This dataframe is then used to
# performs analysis for the particular run or between another runs.


class BenchmarkRun(ReportingObject):
    '''A object to contain benchmark run data which can be used to analyze that run or else to cross-analyze the run with another run.'''

    # Class variables
    amino_acid_details = {}
    CAA, PAA, HAA = set(), set(), set()


    # Human-readable descriptions for the volume breakdown
    by_volume_descriptions = dict(
        SL = 'small-to-large mutations',
        LS = 'large-to-small mutations',
        XX = 'no change in volume',
    )


    def __init__(self, benchmark_run_name, benchmark_run_directory, analysis_directory, dataset_cases, analysis_data, use_single_reported_value,
                 description = None, dataset_description = None, credit = None, take_lowest = 3, generate_plots = True, report_analysis = True, include_derived_mutations = False, recreate_graphs = False, silent = False, burial_cutoff = 0.25,
                 stability_classication_x_cutoff = 1.0, stability_classication_y_cutoff = 1.0, use_existing_benchmark_data = False):
        self.amino_acid_details, self.CAA, self.PAA, self.HAA = BenchmarkRun.get_amino_acid_details()
        self.benchmark_run_name = benchmark_run_name
        self.benchmark_run_directory = benchmark_run_directory
        self.dataset_cases = dataset_cases
        self.analysis_data = analysis_data
        self.analysis_directory = analysis_directory
        self.subplot_directory = os.path.join(self.analysis_directory, self.benchmark_run_name + '_subplots')
        self.analysis_file_prefix = os.path.join(self.analysis_directory, self.benchmark_run_name + '_subplots', self.benchmark_run_name + '_')
        self.use_single_reported_value = use_single_reported_value
        self.description = description
        self.dataset_description = dataset_description
        self.credit = credit
        self.generate_plots = generate_plots
        self.report_analysis = report_analysis
        self.silent = silent
        self.take_lowest = take_lowest
        self.include_derived_mutations = include_derived_mutations
        self.burial_cutoff = burial_cutoff
        self.recreate_graphs = recreate_graphs
        self.stability_classication_x_cutoff = stability_classication_x_cutoff
        self.stability_classication_y_cutoff = stability_classication_y_cutoff
        self.scalar_adjustment = None
        self.analysis_csv_input_filepath = os.path.join(self.benchmark_run_directory, 'analysis_input.csv')
        self.analysis_json_input_filepath = os.path.join(self.benchmark_run_directory, 'analysis_input.json')
        self.analysis_raw_data_input_filepath = os.path.join(self.benchmark_run_directory, 'benchmark_data.json')
        self.analysis_pandas_input_filepath = os.path.join(self.benchmark_run_directory, 'analysis_input.pandas')
        self.metrics_filepath = os.path.join(self.analysis_directory, '{0}_metrics.txt'.format(self.benchmark_run_name))
        self.use_existing_benchmark_data = use_existing_benchmark_data
        self.ddg_analysis_type_description = None
        assert(os.path.exists(self.analysis_csv_input_filepath))
        assert(os.path.exists(self.analysis_json_input_filepath))
        assert(os.path.exists(self.analysis_raw_data_input_filepath))
        if self.generate_plots:
            self.create_subplot_directory() # Create a directory for plots


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


    def report(self, str, fn = None):
        if (not self.silent) and (self.report_analysis):
            if fn:
                fn(str)
            else:
                print(str)


    def create_subplot_directory(self):
        # Create subplot directory
        if not(os.path.exists(self.subplot_directory)):
            try:
                os.mkdir(self.subplot_directory)
                assert(os.path.exists(self.subplot_directory))
            except Exception, e:
                raise colortext.Exception('An exception occurred creating the subplot directory %s.' % self.subplot_directory)


    def create_dataframe(self):
        '''This function creates a dataframe (a matrix with one row per dataset record and one column for fields of interest)
           from the benchmark run and the dataset data.
           For rows with multiple mutations, there may be multiple values for some fields e.g. wildtype residue exposure.
           We take the approach of marking these records as None (to be read as: N/A).
           Another approach is to take averages of continuous and binary values.
           This function also determines a scalar_adjustment used to scale the predictions to try to improve the fraction
           correct score and the MAE.
        '''

        if self.use_existing_benchmark_data and os.path.exists(self.analysis_pandas_input_filepath):
            store = pandas.HDFStore(self.analysis_pandas_input_filepath)
            self.dataframe = store['dataframe']
            self.scalar_adjustment = store['scalar_adjustment'].to_dict()['scalar_adjustment']
            self.ddg_analysis_type = store['ddg_analysis_type'].to_dict()['ddg_analysis_type']
            self.ddg_analysis_type_description = store['ddg_analysis_type_description'].to_dict()['ddg_analysis_type_description']
            store.close()
            return

        analysis_data = self.analysis_data
        dataset_cases = self.dataset_cases
        amino_acid_details, CAA, PAA, HAA = self.amino_acid_details, self.CAA, self.PAA, self.HAA
        burial_cutoff, stability_classication_x_cutoff, stability_classication_y_cutoff = self.burial_cutoff, self.stability_classication_x_cutoff, self.stability_classication_y_cutoff\

        # Create XY data
        self.log('Creating the analysis input file %s and human-readable CSV and JSON versions %s and %s.' % (self.analysis_pandas_input_filepath, self.analysis_csv_input_filepath, self.analysis_json_input_filepath))
        if len(analysis_data) > len(dataset_cases):
            raise colortext.Exception('ERROR: There seems to be an error - there are more predictions than cases in the dataset. Exiting.')
        elif len(analysis_data) < len(dataset_cases):
            self.log('\nWARNING: %d cases missing for analysis; there are %d predictions in the output directory but %d cases in the dataset. The analysis below does not cover the complete dataset.\n' % (len(dataset_cases) - len(analysis_data), len(analysis_data), len(dataset_cases)), colortext.error)

        # ddg_analysis_type can be set to 'DDG' or 'DDG_Top[x]' (e.g. 'DDG_Top3').
        # 'DDG' uses the value reported by ddg_monomer.
        # 'DDG_Top3' (generated by default) uses the metric from Kellogg et al. based on the three lowest scoring mutant structures and the three lowest scoring wildtype structures
        self.ddg_analysis_type = 'DDG_Top%d' % self.take_lowest
        if arguments['--use_single_reported_value']:
            self.ddg_analysis_type = 'DDG'
            self.ddg_analysis_type_description = '\nThe predicted DDG value per case is the single DDG value reported by the ddg_monomer application.'
        elif self.take_lowest == 3:
            self.ddg_analysis_type_description = '\nThe predicted DDG value per case is computed using the {0} lowest-scoring mutant structures and the {0} lowest-scoring wildtype structures as in the paper by Kellogg et al.'.format(self.take_lowest)
        else:
            self.ddg_analysis_type_description = '\nThe predicted DDG value per case is computed using the {0} lowest-scoring mutant structures and the {0} lowest-scoring wildtype structures.'.format(self.take_lowest)
        self.log(self.ddg_analysis_type_description)

        # Initialize the data structures
        csv_headers=[
            'DatasetID', 'PDBFileID', 'Mutations', 'NumberOfMutations', 'Experimental', 'Predicted', 'AbsoluteError', 'StabilityClassification',
            'ResidueCharges', 'VolumeChange',
            'WildTypeDSSPType', 'WildTypeDSSPSimpleSSType', 'WildTypeDSSPExposure',
            'WildTypeSCOPClass', 'WildTypeSCOPFold', 'WildTypeSCOPClassification',
            'WildTypeExposure', 'WildTypeAA', 'MutantAA', 'HasGPMutation',
            'PDBResolution', 'PDBResolutionBin', 'MonomerLength', 'NumberOfDerivativeErrors',
        ]
        csv_file = [','.join(csv_headers)]

        # Set the PDB input path
        pdb_data = {}
        try:
            pdb_data_ = json.loads(read_file('../../input/json/pdbs.json'))
            for k, v in pdb_data_.iteritems():
                pdb_data[k.upper()] = v
        except Exception, e:
            self.log('input/json/pdbs.json could not be found - PDB-specific analysis cannot be performed.', colortext.error)

        # Create the dataframe
        for record_id, predicted_data in sorted(analysis_data.iteritems()):

            record = dataset_cases[record_id]

            # Ignore derived mutations if appropriate
            if record['DerivedMutation'] and not self.include_derived_mutations:
                continue

            # Initialize variables. For ambiguous cases where the set of distinct values has multiple values, we default to None
            residue_charge, residue_charges = None, set()
            exposure, exposures = None, set()
            volume_change, volume_changes = None, set()
            record_wtaa, wtaas = None, set()
            record_mutaa, mutaas = None, set()
            DSSPSimpleSSType, DSSPSimpleSSTypes = None, set()
            DSSPType, DSSPTypes = None, set()
            DSSPExposure, DSSPExposures = None, set()
            scops = set()
            pdb_chains = set()
            mutation_string = []
            num_derivative_errors = predicted_data.get('Errors', {}).get('Derivative error count', 0)

            mutations = record['Mutations']
            for m in mutations:

                wtaa = m['WildTypeAA']
                mutaa = m['MutantAA']
                mutation_string.append('{0} {1}{2}{3}'.format(m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA']))

                # Residue types and chain
                wtaas.add(wtaa)
                mutaas.add(mutaa)
                pdb_chains.add(m['Chain'])
                scops.add(m['SCOP class'])
                DSSPSimpleSSTypes.add(m['DSSPSimpleSSType'])
                DSSPTypes.add(m['DSSPType'])
                DSSPExposures.add(m['DSSPExposure'])

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

            # Create a string representing the mutations (useful for labeling rather than analysis)
            mutation_string = '; '.join(mutation_string)

            # Taking unique values, determine the residue charges of the wildtype and mutant residues, the wildtype residue exposure, and the relative change in van der Waals volume
            if len(residue_charges) == 1: residue_charge = residue_charges.pop()
            if len(exposures) == 1: exposure = exposures.pop()
            if len(volume_changes) == 1: volume_change = volume_changes.pop()

            # Taking unique values, determine the wildtype and mutant residue types
            all_residues = wtaas.union(mutaas)
            if len(wtaas) == 1: record_wtaa = wtaas.pop()
            if len(mutaas) == 1: record_mutaa = mutaas.pop()
            assert(len(pdb_chains) == 1) # we expect monomeric cases
            pdb_chain = pdb_chains.pop()

            # Taking unique values, determine the secondary structure and residue exposures from the DSSP data in the dataset
            if len(DSSPSimpleSSTypes) == 1: DSSPSimpleSSType = DSSPSimpleSSTypes.pop()
            if len(DSSPTypes) == 1: DSSPType = DSSPTypes.pop()
            if len(DSSPExposures) == 1: DSSPExposure = DSSPExposures.pop()

            # Determine the SCOP classification from the SCOPe data in the dataset
            full_scop_classification, scop_class, scop_fold = None, None, None
            if len(scops) > 1:
                self.log('Warning: There is more than one SCOPe class for record {0}.'.format(record_id), colortext.warning)
            else:
                full_scop_classification = scops.pop()
                scop_tokens = full_scop_classification.split('.')
                scop_class = scop_tokens[0]
                if len(scop_tokens) > 1:
                    scop_fold = '.'.join(scop_tokens[0:2])

            # Calculate the stability classification and absolute_error for this case
            stability_classification = fraction_correct([record['DDG']], [predicted_data[self.ddg_analysis_type]], x_cutoff = stability_classication_x_cutoff, y_cutoff = stability_classication_y_cutoff)
            assert(stability_classification == 0 or stability_classification == 1)
            #stability_classification = int(stability_classification)
            absolute_error = abs(record['DDG'] - predicted_data[self.ddg_analysis_type])

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

            # Mark mutations involving glycine or proline
            has_gp_mutation = 'G' in all_residues or 'P' in all_residues

            # Create the data matrix
            dataframe_record = dict(
                DatasetID = record_id,
                PDBFileID = record['PDBFileID'],
                Mutations = mutation_string,
                NumberOfMutations = len(mutations),
                Experimental = record['DDG'],
                Predicted = predicted_data[self.ddg_analysis_type],
                AbsoluteError = absolute_error,
                StabilityClassification = stability_classification,
                ResidueCharges = residue_charge,
                VolumeChange = volume_change,
                HasGPMutation = int(has_gp_mutation),
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
                NumberOfDerivativeErrors = num_derivative_errors,
                )

            for h in csv_headers:
                assert(',' not in str(dataframe_record[h]))
            csv_file.append(','.join([str(dataframe_record[h]) for h in csv_headers]))
            assert(sorted(csv_headers) == sorted(dataframe_record.keys()))

        # Create the CSV file in memory (we are not done added data just yet) and pass it to pandas
        dataframe = pandas.read_csv(StringIO.StringIO('\n'.join(csv_file)), sep=',', header=0, skip_blank_lines=True, index_col = 0)
        self.dataframe = dataframe

        # Report the SCOPe classification counts
        SCOP_classifications = set(dataframe['WildTypeSCOPClassification'].values.tolist())
        SCOP_folds = set(dataframe['WildTypeSCOPFold'].values.tolist())
        SCOP_classes = set(dataframe['WildTypeSCOPClass'].values.tolist())
        self.log('The mutated residues span {0} unique SCOP(e) classifications in {1} unique SCOP(e) folds and {2} unique SCOP(e) classes.'.format(len(SCOP_classifications), len(SCOP_folds), len(SCOP_classes)), colortext.message)

        # Plot the optimum y-cutoff over a range of x-cutoffs for the fraction correct metric. Include the user's cutoff in the range
        self.log('Determining a scalar adjustment with which to scale the predicted values to improve the fraction correct measurement.', colortext.warning)
        self.scalar_adjustment, plot_filename = self.plot_optimum_prediction_fraction_correct_cutoffs_over_range(min(self.stability_classication_x_cutoff, 0.5), max(self.stability_classication_x_cutoff, 3.0), suppress_plot = True)

        # Add new columns derived from the adjusted values
        dataframe['Predicted_adj'] = dataframe['Predicted'] / self.scalar_adjustment
        dataframe['AbsoluteError_adj'] = (dataframe['Experimental'] - dataframe['Predicted_adj']).abs()
        add_fraction_correct_values_to_dataframe(dataframe, 'Experimental', 'Predicted_adj', 'StabilityClassification_adj',  x_cutoff = stability_classication_x_cutoff, y_cutoff = stability_classication_y_cutoff)

        # Write the dataframe out to CSV
        dataframe.to_csv(self.analysis_csv_input_filepath, sep = ',', header = True)

        # Write the dataframe out to JSON
        # Note: I rolled my own as dataframe.to_dict(orient = 'records') gives us the correct format but discards the DatasetID (index) field
        json_records = {}
        indices = dataframe.index.values.tolist()
        for i in indices:
            json_records[i] = {}
        for k, v in dataframe.to_dict().iteritems():
            for i, v in v.iteritems():
                assert(k not in json_records[i])
                json_records[i][k] = v
        write_file(self.analysis_json_input_filepath, json.dumps(json_records, indent = 4, sort_keys=True))

        # Write the values computed in this function out to disk
        if os.path.exists(self.analysis_pandas_input_filepath):
            os.remove(self.analysis_pandas_input_filepath)
        store = pandas.HDFStore(self.analysis_pandas_input_filepath)
        store['dataframe'] = dataframe
        store['scalar_adjustment'] = pandas.Series(dict(scalar_adjustment = self.scalar_adjustment))
        store['ddg_analysis_type'] = pandas.Series(dict(ddg_analysis_type = self.ddg_analysis_type))
        store['ddg_analysis_type_description'] = pandas.Series(dict(ddg_analysis_type_description = self.ddg_analysis_type_description))
        store.close()


    def analyze(self):
        '''This function runs the analysis and creates the plots and summary file.'''
        self.calculate_metrics()
        self.plot()


    def calculate_metrics(self):
        '''Calculates the main metrics for the benchmark run and writes them to file.'''

        dataframe = self.dataframe
        scalar_adjustment = self.scalar_adjustment

        metrics_textfile = [self.ddg_analysis_type_description]

        if self.include_derived_mutations:
            metrics_textfile.append('\nRunning analysis (derived mutations will be included):')
        else:
            metrics_textfile.append('\nRunning analysis (derived mutations will be omitted):')
        self.report(metrics_textfile[-1], fn = colortext.message)
        metrics_textfile.append('The stability classification cutoffs are: Experimental=%0.2f kcal/mol, Predicted=%0.2f energy units.' % (self.stability_classication_x_cutoff, self.stability_classication_y_cutoff))
        self.report(metrics_textfile[-1], fn = colortext.warning)

        amino_acid_details, CAA, PAA, HAA = self.amino_acid_details, self.CAA, self.PAA, self.HAA

        # This dict is used for the print-statement below
        volume_groups = {}
        for aa_code, aa_details in amino_acid_details.iteritems():
            v = int(aa_details['van der Waals volume']) # Note: I only convert to int here to match the old script behavior and because all volumes are integer values so it does not do any harm
            volume_groups[v] = volume_groups.get(v, [])
            volume_groups[v].append(aa_code)

        metrics_textfile.append('\n\nSection 1. Breakdown by volume.')
        metrics_textfile.append('A case is considered a small-to-large (resp. large-to-small) mutation if all of the wildtype residues have a smaller (resp. larger) van der Waals volume than the corresponding mutant residue. The order is defined as %s so some cases are considered to have no change in volume e.g. MET -> LEU.' % (' < '.join([''.join(sorted(v)) for k, v in sorted(volume_groups.iteritems())])))
        self.report('\n'.join(metrics_textfile[-2:]), fn = colortext.sprint)
        for subcase in ('XX', 'SL', 'LS'):
            subcase_dataframe = dataframe[dataframe['VolumeChange'] == subcase]
            metrics_textfile.append('\n' + '*'*10 + (' Statistics - %s (%d cases)' % (BenchmarkRun.by_volume_descriptions[subcase], len(subcase_dataframe))) +'*'*10)
            metrics_textfile.append(format_stats_for_printing(get_xy_dataset_statistics_pandas(subcase_dataframe, 'Experimental', 'Predicted', fcorrect_x_cutoff = self.stability_classication_x_cutoff, fcorrect_y_cutoff = self.stability_classication_y_cutoff)))
            self.report('\n'.join(metrics_textfile[-2:]), fn = colortext.sprint)

        metrics_textfile.append('\n\nSection 2. Separating out mutations involving glycine or proline.')
        metrics_textfile.append('This cases may involve changes to secondary structure so we separate them out here.')
        subcase_dataframe = dataframe[dataframe['HasGPMutation'] == 1]
        metrics_textfile.append('\n' + '*'*10 + (' Statistics - cases with G or P (%d cases)' % (len(subcase_dataframe))) +'*'*10)
        metrics_textfile.append(format_stats_for_printing(get_xy_dataset_statistics_pandas(subcase_dataframe, 'Experimental', 'Predicted', fcorrect_x_cutoff = self.stability_classication_x_cutoff, fcorrect_y_cutoff = self.stability_classication_y_cutoff)))
        self.report('\n'.join(metrics_textfile[-4:]), fn = colortext.sprint)
        subcase_dataframe = dataframe[dataframe['HasGPMutation'] == 0]
        metrics_textfile.append('\n' + '*'*10 + (' Statistics - cases without G or P (%d cases)' % (len(subcase_dataframe))) +'*'*10)
        metrics_textfile.append(format_stats_for_printing(get_xy_dataset_statistics_pandas(subcase_dataframe, 'Experimental', 'Predicted', fcorrect_x_cutoff = self.stability_classication_x_cutoff, fcorrect_y_cutoff = self.stability_classication_y_cutoff)))
        self.report('\n'.join(metrics_textfile[-2:]), fn = colortext.sprint)


        metrics_textfile.append('\nSection 3. Entire dataset using a scaling factor of 1/%.03f to improve the fraction correct metric.' % scalar_adjustment)
        self.report(metrics_textfile[-1], fn = colortext.sprint)
        metrics_textfile.append('Warning: Results in this section use an averaged scaling factor to improve the value for the fraction correct metric. This scalar will vary over benchmark runs so these results should not be interpreted as performance results; they should be considered as what could be obtained if the predicted values were scaled by a "magic" value.')
        self.report(metrics_textfile[-1], fn = colortext.warning)
        metrics_textfile.append('\n' + '*'*10 + (' Statistics - complete dataset (%d cases)' % len(dataframe)) +'*'*10)
        # For these statistics, we assume that we have reduced any scaling issues and use the same cutoff for the Y-axis as the user specified for the X-axis
        metrics_textfile.append(format_stats_for_printing(get_xy_dataset_statistics_pandas(dataframe, 'Experimental', 'Predicted_adj', fcorrect_x_cutoff = self.stability_classication_x_cutoff, fcorrect_y_cutoff = self.stability_classication_x_cutoff)))
        self.report('\n'.join(metrics_textfile[-2:]), fn = colortext.sprint)

        metrics_textfile.append('\n\nSection 4. Entire dataset.')
        self.report(metrics_textfile[-1], fn = colortext.sprint)
        metrics_textfile.append('Overall statistics')
        self.report(metrics_textfile[-1], fn = colortext.message)
        metrics_textfile.append('\n' + '*'*10 + (' Statistics - complete dataset (%d cases)' % len(dataframe)) +'*'*10)
        metrics_textfile.append(format_stats_for_printing(get_xy_dataset_statistics_pandas(dataframe, 'Experimental', 'Predicted', fcorrect_x_cutoff = self.stability_classication_x_cutoff, fcorrect_y_cutoff = self.stability_classication_y_cutoff)))
        self.report('\n'.join(metrics_textfile[-2:]), fn = colortext.sprint)

        # There is probably a better way of writing the pandas code here
        record_with_most_errors = (dataframe[['PDBFileID', 'NumberOfDerivativeErrors', 'Mutations']].sort('NumberOfDerivativeErrors')).tail(1)
        record_index = record_with_most_errors.index.tolist()[0]
        pdb_id, num_errors, mutation_str = dataframe.loc[record_index, 'PDBFileID'], dataframe.loc[record_index, 'NumberOfDerivativeErrors'], dataframe.loc[record_index, 'Mutations']
        if num_errors > 0:
            metrics_textfile.append('\n\nDerivative errors were found in the run. Record #{0} - {1}, {2} - has the most amount ({3}) of derivative errors.'.format(record_index, pdb_id, mutation_str, num_errors))
            self.report(metrics_textfile[-1], fn = colortext.warning)

        # Write the analysis to file
        write_file(self.metrics_filepath, '\n'.join(metrics_textfile))


    def plot(self):

        if not self.generate_plots:
            return

        dataframe = self.dataframe
        graph_order = []

        # Create a subtitle for the first page
        subtitle = self.benchmark_run_name
        if self.description:
            tokens = self.description.split(' ')
            description_lines = []
            s = []
            for t in tokens:
                if s and (len(s) + len(t) + 1 > 32):
                    description_lines.append(' '.join(s))
                    s = [t]
                else:
                    s.append(t)
            if s: description_lines.append(' '.join(s))
            subtitle += '\n{0}'.format('\n'.join(description_lines))
        if self.dataset_description:
            subtitle += '\n{0}'.format(self.dataset_description)

        # Plot which y-cutoff yields the best value for the fraction correct metric
        scalar_adjustment, scalar_adjustment_calculation_plot = self.plot_optimum_prediction_fraction_correct_cutoffs_over_range(min(self.stability_classication_x_cutoff, 0.5), max(self.stability_classication_x_cutoff, 3.0), suppress_plot = False)
        assert(self.scalar_adjustment == scalar_adjustment)

        # Plot which the optimum y-cutoff given the specified or default x-cutoff
        optimal_predictive_cutoff_plot = self.plot_optimum_prediction_fraction_correct_cutoffs(self.stability_classication_x_cutoff)

        # Create a scatterplot for the results
        main_scatterplot = '{0}main_scatterplot.png'.format(self.analysis_file_prefix)
        if not(os.path.exists(main_scatterplot) and not(self.recreate_graphs)):
            self.log('Saving scatterplot to %s.' % main_scatterplot)
            plot_pandas(dataframe, 'Experimental', 'Predicted', main_scatterplot, RInterface.correlation_coefficient_gplot, title = 'Experimental vs. Prediction')

        graph_order.append(self.create_section_slide('{0}section_1.png'.format(self.analysis_file_prefix), 'Main metrics', subtitle, self.credit))
        graph_order.append(main_scatterplot)

        # Plot a histogram of the absolute errors
        graph_order.append(self.plot_absolute_error_histogram('{0}absolute_errors'.format(self.analysis_file_prefix), 'AbsoluteError'))
        graph_order.append(self.create_section_slide('{0}section_2.png'.format(self.analysis_file_prefix), 'Adjustments', 'Optimization of the cutoffs\nfor the fraction correct metric'))
        graph_order.append(scalar_adjustment_calculation_plot)
        graph_order.append(optimal_predictive_cutoff_plot)

        # Create a scatterplot and histogram for the adjusted results
        main_adj_scatterplot = '{0}main_adjusted_with_scalar_scatterplot.png'.format(self.analysis_file_prefix)
        if not(os.path.exists(main_adj_scatterplot) and not(self.recreate_graphs)):
            self.log('Saving scatterplot to %s.' % main_adj_scatterplot)
            plot_pandas(dataframe, 'Experimental', 'Predicted_adj', main_adj_scatterplot, RInterface.correlation_coefficient_gplot, title = 'Experimental vs. Prediction: adjusted scale')
        graph_order.append(main_adj_scatterplot)
        graph_order.append(self.plot_absolute_error_histogram('{0}absolute_errors_adjusted_with_scalar'.format(self.analysis_file_prefix), 'AbsoluteError_adj'))

        # Scatterplots colored by residue context / change on mutation
        graph_order.append(self.create_section_slide('{0}section_3.png'.format(self.analysis_file_prefix), 'Residue context'))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Residue charges', self.scatterplot_charges, '{0}scatterplot_charges.png'.format(self.analysis_file_prefix)))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Exposure (cutoff = %0.2f)' % self.burial_cutoff, self.scatterplot_exposure, '{0}scatterplot_exposure.png'.format(self.analysis_file_prefix)))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Change in volume', self.scatterplot_volume, '{0}scatterplot_volume.png'.format(self.analysis_file_prefix)))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Wildtype residue s.s.', self.scatterplot_ss, '{0}scatterplot_ss.png'.format(self.analysis_file_prefix)))

        # Scatterplots colored by SCOPe classification
        graph_order.append(self.create_section_slide('{0}section_4.png'.format(self.analysis_file_prefix), 'SCOPe classes'))
        SCOP_classifications = set(dataframe['WildTypeSCOPClassification'].values.tolist())
        SCOP_folds = set(dataframe['WildTypeSCOPFold'].values.tolist())
        SCOP_classes = set(dataframe['WildTypeSCOPClass'].values.tolist())
        if len(SCOP_classes) <= 25:
            graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - WT residue SCOP class', self.scatterplot_scop_class, '{0}scatterplot_scop_class.png'.format(self.analysis_file_prefix)))
        if len(SCOP_folds) <= 25:
            graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - WT residue SCOP fold', self.scatterplot_scop_fold, '{0}scatterplot_scop_fold.png'.format(self.analysis_file_prefix)))
        if len(SCOP_classifications) <= 25:
            graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - WT residue SCOP classification', self.scatterplot_scop_classification, '{0}scatterplot_scop_classification.png'.format(self.analysis_file_prefix)))

        # Scatterplots colored by residue types
        graph_order.append(self.create_section_slide('{0}section_5.png'.format(self.analysis_file_prefix), 'Residue types'))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Wildtype', self.scatterplot_wildtype_aa, '{0}scatterplot_wildtype_aa.png'.format(self.analysis_file_prefix)))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Mutant', self.scatterplot_mutant_aa, '{0}scatterplot_mutant_aa.png'.format(self.analysis_file_prefix)))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Glycine/Proline', self.scatterplot_GP, '{0}scatterplot_gp.png'.format(self.analysis_file_prefix)))

        # Scatterplots colored PDB resolution and chain length
        graph_order.append(self.create_section_slide('{0}section_6.png'.format(self.analysis_file_prefix), 'Chain properties'))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - PDB resolution', self.scatterplot_pdb_res_binned, '{0}scatterplot_pdb_res_binned.png'.format(self.analysis_file_prefix)))
        graph_order.append(self.scatterplot_generic('Experimental vs. Prediction - Chain length', self.scatterplot_chain_length, '{0}scatterplot_chain_length.png'.format(self.analysis_file_prefix)))

        # Errors / debugging
        graph_order.append(self.create_section_slide('{0}section_7.png'.format(self.analysis_file_prefix), 'Errors / debugging'))
        graph_order.append(self.plot_derivative_error_barchart())

        # Make sure all of the graphs have been created
        relative_graph_paths = [os.path.join(self.benchmark_run_name + '_subplots', os.path.split(g)[1]) for g in graph_order]
        for rp in relative_graph_paths:
            assert(os.path.exists(os.path.join(self.analysis_directory, rp)))

        # Copy the analysis input files into the analysis directory - these files are duplicated but it makes it easier to share data
        shutil.copyfile(self.analysis_csv_input_filepath, os.path.join(self.analysis_directory, self.benchmark_run_name + '_' + os.path.split(self.analysis_csv_input_filepath)[1]))
        shutil.copyfile(self.analysis_json_input_filepath, os.path.join(self.analysis_directory, self.benchmark_run_name + '_' + os.path.split(self.analysis_json_input_filepath)[1]))
        shutil.copyfile(self.analysis_raw_data_input_filepath, os.path.join(self.analysis_directory, self.benchmark_run_name + '_' + os.path.split(self.analysis_raw_data_input_filepath)[1]))

        # Combine the plots into a PDF file
        plot_file = os.path.join(self.analysis_directory, '{0}_benchmark_plots.pdf'.format(self.benchmark_run_name))
        self.log('\n\nCreating a PDF containing all of the plots: {0}'.format(plot_file), colortext.message)
        try:
            if os.path.exists(plot_file):
                # raise Exception('remove this exception') #todo
                os.remove(plot_file)
            p = subprocess.Popen(shlex.split('convert {0} {1}'.format(' '.join(relative_graph_paths), plot_file)), cwd = self.analysis_directory)
            stdoutdata, stderrdata = p.communicate()
            if p.returncode != 0: raise Exception('')
        except:
            self.log('An error occurred while combining the positional scatterplots using the convert application (ImageMagick).', colortext.error)


    def compare(self, other):
        '''Compare this benchmark run with another run.'''
        pass
        #todo implement this


    def determine_optimum_fraction_correct_cutoffs(self, dataframe, stability_classication_x_cutoff):
        '''Determines the value of stability_classication_y_cutoff which approximately maximizes the fraction correct
           measurement w.r.t. a fixed stability_classication_x_cutoff. This function uses discrete sampling and so it
           may miss the actual maximum. We use two rounds of sampling: i) a coarse-grained sampling (0.1 energy unit
           intervals); and ii) finer sampling (0.01 unit intervals).
           In both rounds, we choose the one corresponding to a lower value for the cutoff in cases of multiple maxima.'''

        # Determine the value for the fraction correct y-value (predicted) cutoff which will approximately yield the
        # maximum fraction-correct value

        fraction_correct_range = []

        # Round 1 : Coarse sampling. Test 0.5 -> 8.0 in 0.1 increments
        for z in range(5, 80):
            w = float(z) / 10.0
            fraction_correct_range.append((w, fraction_correct_pandas(dataframe, 'Experimental', 'Predicted', x_cutoff = stability_classication_x_cutoff, y_cutoff = w)))

        max_value_cutoff, max_value = fraction_correct_range[0][0], fraction_correct_range[0][1]
        for p in fraction_correct_range:
            if p[1] > max_value:
                max_value_cutoff, max_value = p[0], p[1]

        # Round 2 : Finer sampling. Test max_value_cutoff - 0.1 -> max_value_cutoff + 0.1 in 0.01 increments
        for z in range(int((max_value_cutoff - 0.1) * 100), int((max_value_cutoff + 0.1) * 100)):
            w = float(z) / 100.0
            fraction_correct_range.append((w, fraction_correct_pandas(dataframe, 'Experimental', 'Predicted', x_cutoff = stability_classication_x_cutoff, y_cutoff = w)))
        fraction_correct_range = sorted(set(fraction_correct_range)) # sort so that we find the lowest cutoff value in case of duplicate fraction correct values
        max_value_cutoff, max_value = fraction_correct_range[0][0], fraction_correct_range[0][1]
        for p in fraction_correct_range:
            if p[1] > max_value:
                max_value_cutoff, max_value = p[0], p[1]

        return max_value_cutoff, max_value, fraction_correct_range


    def create_section_slide(self, plot_filename, title, subtitle = '', footer = ''):

        if os.path.exists(plot_filename) and not(self.recreate_graphs):
            return plot_filename

        R_filename = os.path.splitext(plot_filename)[0] + '.R'

        title_size = 8
        longest_line = max(map(len, title.split('\n')))
        if longest_line <= 11:
            title_size = 16
        elif longest_line <= 16:
            title_size = 12
        elif longest_line <= 19:
            title_size = 10

        subtitle_size = 6
        longest_line = max(map(len, subtitle.split('\n')))
        if longest_line <= 11:
            subtitle_size = 16
        elif longest_line <= 16:
            subtitle_size = 12
        elif longest_line <= 19:
            subtitle_size = 10
        elif longest_line <= 32:
            subtitle_size = 8
        subtitle_size = min(title_size - 2, subtitle_size)

        footer_size = 6

        # Create plot
        if self.generate_plots:
            r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)

p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() +
     coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
     theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) +
     geom_text(hjust=0, vjust=1, size=%(title_size)d, color="black", aes(5, 90, fontface="plain", family = "Times", face = "bold", label="%(title)s")) +
     geom_text(hjust=0, vjust=1, size=%(subtitle_size)d, color="black", aes(5, 50, fontface="plain", family = "Times", face = "bold", label="%(subtitle)s")) +
     geom_text(hjust=0, vjust=1, size=%(footer_size)d, color="black", aes(5, 5, fontface="italic", family = "Times", face = "bold", label="%(footer)s"))
p
dev.off()'''
            self.log('Section slide %s.' % plot_filename)
            RInterface._runRScript(r_script % locals())
            return plot_filename


    def plot_optimum_prediction_fraction_correct_cutoffs(self, stability_classication_x_cutoff):

        # Determine the optimal values
        max_value_cutoff, max_value, fraction_correct_range = self.determine_optimum_fraction_correct_cutoffs(self.dataframe, stability_classication_x_cutoff)

        # Filenames
        output_filename_prefix = '{0}optimum_fraction_correct_at_{1}_kcal_mol'.format(self.analysis_file_prefix, '%.2f' % stability_classication_x_cutoff)
        plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'

        if os.path.exists(plot_filename) and not(self.recreate_graphs):
            return plot_filename

        # Create CSV input
        lines = ['NeutralityCutoff,FractionCorrect,C']
        for p in fraction_correct_range:
            if p[1] == max_value:
                lines.append(','.join(map(str, (p[0], p[1], 'best'))))
            else:
                lines.append(','.join(map(str, (p[0], p[1], 'other'))))
        write_file(csv_filename, '\n'.join(lines))

        # Create plot
        if self.generate_plots:
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
            self.log('Saving plot of approximate optimal fraction correct cutoffs to %s.' % plot_filename)
            RInterface._runRScript(r_script % locals())
            return plot_filename


    def plot_optimum_prediction_fraction_correct_cutoffs_over_range(self, min_stability_classication_x_cutoff, max_stability_classication_x_cutoff, suppress_plot = False):
        '''Plots the optimum cutoff for the predictions to maximize the fraction correct metric over a range of experimental cutoffs.
           Returns the average scalar corresponding to the best value of fraction correct over a range of cutoff values for the experimental cutoffs.'''

        # Filenames
        output_filename_prefix = '{0}optimum_fraction_correct_at_varying_kcal_mol'.format(self.analysis_file_prefix)
        plot_filename = None
        if not suppress_plot:
            plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'

        # Create CSV input
        lines = ['ExperimentalCutoff,BestPredictionCutoff']
        x_cutoff = min_stability_classication_x_cutoff
        x_values = []
        y_values = []
        avg_scale = 0
        plot_graph = self.generate_plots and not(suppress_plot)
        while x_cutoff < max_stability_classication_x_cutoff + 0.1:
            max_value_cutoff, max_value, fraction_correct_range = self.determine_optimum_fraction_correct_cutoffs(self.dataframe, x_cutoff)
            if plot_graph:
                lines.append(','.join(map(str, (x_cutoff, max_value_cutoff))))
            x_values.append(x_cutoff)
            y_values.append(max_value_cutoff)
            avg_scale += max_value_cutoff / x_cutoff
            x_cutoff += 0.1

        if plot_graph:
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
        if plot_graph:
            if not(os.path.exists(plot_filename) and not(self.recreate_graphs)):
                self.log('Saving scatterplot to %s.' % plot_filename)
                self.log('Saving plot of approximate optimal fraction correct cutoffs over varying experimental cutoffs to %s.' % plot_filename)

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

        return average_scalar, plot_filename


    def plot_derivative_error_barchart(self):

        # Filenames
        output_filename_prefix = '{0}errors_by_pdb_id'.format(self.analysis_file_prefix)
        plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'

        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['PDBFileID', 'NumberOfDerivativeErrors']]
        new_dataframe.columns = ['PDB', 'AverageDerivativeErrorCount']
        new_dataframe = new_dataframe.groupby(['PDB'])['AverageDerivativeErrorCount'].mean()
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        if os.path.exists(plot_filename) and not(self.recreate_graphs):
            return plot_filename

        # Create plot
        firebrick = plot_colors['firebrick']
        brown = plot_colors['brown']
        self.log('Saving barchart to %s.' % plot_filename)
        title = 'Average count of Inaccurate G! errors by PDB ID'
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

b <- ggplot(plot_data, aes(x=PDB, y=AverageDerivativeErrorCount)) +
     geom_bar(stat='identity', colour = "%(brown)s", fill = "%(firebrick)s") +
     ggtitle("%(title)s") +
     xlab("PDB ID") +
     ylab("Derivative errors (average)") +
     coord_flip()
b

#m <- ggplot(plot_data, aes(x=AbsoluteError)) +
#    geom_histogram(colour = "%(brown)s", fill = "%(firebrick)s", binwidth = 0.5) +
#    ggtitle("%(title)s") +
#    xlab("Absolute error (kcal/mol - energy units)") +
#    ylab("Number of cases")
#m
dev.off()'''
        RInterface._runRScript(r_script % locals())
        return plot_filename


    def plot_absolute_error_histogram(self, output_filename_prefix, data_series):

        # Filenames
        plot_filename = output_filename_prefix + '.png'
        csv_filename = output_filename_prefix + '.txt'
        R_filename = output_filename_prefix + '.R'

        if os.path.exists(plot_filename) and not(self.recreate_graphs):
            return plot_filename

        # Create CSV input
        new_dataframe = self.dataframe[[data_series]]
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        if not self.generate_plots:
            return

        # Create plot
        self.log('Saving scatterplot to %s.' % plot_filename)
        title = 'Distribution of absolute errors (prediction - observed)'
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

m <- ggplot(plot_data, aes(x=%(data_series)s)) +
    geom_histogram(colour = "darkgreen", fill = "green", binwidth = 0.5) +
    ggtitle("%(title)s") +
    xlab("Absolute error (kcal/mol - energy units)") +
    ylab("Number of cases")
m
dev.off()'''
        RInterface._runRScript(r_script % locals())
        return plot_filename


    def scatterplot_generic(self, title, plotfn, plot_filename):
        if os.path.exists(plot_filename) and not(self.recreate_graphs):
            return plot_filename

        csv_filename = os.path.splitext(plot_filename)[0] + '.txt'
        plot_commands = plotfn(title, csv_filename)
        r_script = '''library(ggplot2)
library(gridExtra)
library(scales)
library(qualV)

png('%(plot_filename)s', height=4096, width=4096, bg="white", res=600)
plot_data <- read.csv('%(csv_filename)s', header=T)

%(plot_commands)s

dev.off()''' % locals()
        if self.generate_plots:
            self.log('Saving scatterplot to %s.' % plot_filename)
            RInterface._runRScript(r_script)
            return plot_filename


    def scatterplot_color_by_series(self, colorseries, xseries = "Experimental", yseries = "Predicted", title = '', plot_scale = '', point_opacity = 0.4, extra_commands = ''):

        # Compute MAE
        mae_str = ''
        if xseries == 'Experimental':
            if yseries == 'Predicted':
                mae_str = self.dataframe['AbsoluteError'].mean()
            elif yseries == 'Predicted_adj':
                mae_str = self.dataframe['AbsoluteError_adj'].mean()
        if mae_str:
            mae_str = 'MAE = {0:.3f}'.format(mae_str)

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
    geom_text(hjust=0, size=4, aes(xpos, ypos_mae, fontface="plain", family = "sans", label="%(mae_str)s"))
p

''' % locals()


    def scatterplot_charges(self, title, csv_filename):
        '''Scatterplot by residue charge.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'ResidueCharges']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    values = c( "None" = '#777777', "Change" = '%(cornflower_blue)s', "Polar/Charged" = 'magenta', "Hydrophobic/Non-polar" = 'green'),
    labels = c( "None" = "N/A", "Change" = "Change", "Polar/Charged" = "Polar/Charged", "Hydrophobic/Non-polar" = "Hydrophobic/Non-polar"))''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "ResidueCharges", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_exposure(self, title, csv_filename):
        '''Scatterplot by exposure class.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'WildTypeExposure']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.columns = ['Experimental', 'Predicted', 'Exposure', 'Opacity'] # rename the exposure column
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    values = c( "None" = '#777777', "B" = '%(brown)s', "E" = '%(purple)s'),
    labels = c( "None" = "N/A", "B" = "Buried", "E" = "Exposed"))''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "Exposure", title = title, plot_scale = plot_scale)


    def scatterplot_volume(self, title, csv_filename):
        '''Scatterplot by change in volume upon mutation.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'VolumeChange']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    values = c( "None" = '#777777', "SL" = '%(brown)s', "LS" = '%(purple)s', 'XX' = "%(cornflower_blue)s"),
    labels = c( "None" = "N/A", "SL" = "Increase", "LS" = "Decrease", "XX" = "No change"))''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "VolumeChange", title = title, plot_scale = plot_scale)


    def scatterplot_ss(self, title, csv_filename):
        '''Scatterplot by secondary structure.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'WildTypeDSSPSimpleSSType']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.columns = ['Experimental', 'Predicted', 'WTSecondaryStructure', 'Opacity'] # rename the s.s. column
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    name="Secondary structure",
    values = c( "None" = '#777777', "H" = 'magenta', "S" = 'orange', "O" = '%(cornflower_blue)s'),
    labels = c( "None" = "N/A", "H" = "Helix", "S" = "Sheet", "O" = "Other"))''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "WTSecondaryStructure", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_scop_class(self, title, csv_filename):
        '''Scatterplot by SCOPe class.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'WildTypeSCOPClass']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        return self.scatterplot_color_by_series(colorseries = "WildTypeSCOPClass", title = title, point_opacity = 0.6)


    def scatterplot_scop_fold(self, title, csv_filename):
        '''Scatterplot by SCOPe fold.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'WildTypeSCOPFold']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        return self.scatterplot_color_by_series(colorseries = "WildTypeSCOPFold", title = title, point_opacity = 0.6)


    def scatterplot_scop_classification(self, title, csv_filename):
        '''Scatterplot by SCOPe classification.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'WildTypeSCOPClassification']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        return self.scatterplot_color_by_series(colorseries = "WildTypeSCOPClassification", title = title, point_opacity = 0.6)


    def scatterplot_wildtype_aa(self, title, csv_filename):
        '''Scatterplot by wildtype residue.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'WildTypeAA']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    name="Residue",
    values = c( "None" = '#808080', "A" = '#FF0000', "C" = '#BFBF00', "D" = '#008000', "E" = "#80FFFF", "F" = "#8080FF", "G" = "#BF40BF", "H" = "#A0A424", "I" = "#411BEA", "K" = "#1EAC41", "L" = "#F0C80E", "M" = "#B430E5", "N" = "#ED7651", "P" = "#19CB97", "Q" = "#362698", "R" = "#7E7EB8", "S" = "#603000", "T" = "#A71818", "V" = "#DF8020", "W" = "#E75858", "Y" = "#082008"),
    labels = c( "None" = "N/A", "A" = "A", "C" = "C", "D" = "D", "E" = "E", "F" = "F", "G" = "G", "H" = "H", "I" = "I", "K" = "K", "L" = "L", "M" = "M", "N" = "N", "P" = "P", "Q" = "Q", "R" = "R", "S" = "S", "T" = "T", "V" = "V", "W" = "W", "Y" = "Y"))
    ''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "WildTypeAA", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_mutant_aa(self, title, csv_filename):
        '''Scatterplot by mutant residue.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'MutantAA']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    name="Residue",
    values = c( "None" = '#808080', "A" = '#FF0000', "C" = '#BFBF00', "D" = '#008000', "E" = "#80FFFF", "F" = "#8080FF", "G" = "#BF40BF", "H" = "#A0A424", "I" = "#411BEA", "K" = "#1EAC41", "L" = "#F0C80E", "M" = "#B430E5", "N" = "#ED7651", "P" = "#19CB97", "Q" = "#362698", "R" = "#7E7EB8", "S" = "#603000", "T" = "#A71818", "V" = "#DF8020", "W" = "#E75858", "Y" = "#082008"),
    labels = c( "None" = "N/A", "A" = "A", "C" = "C", "D" = "D", "E" = "E", "F" = "F", "G" = "G", "H" = "H", "I" = "I", "K" = "K", "L" = "L", "M" = "M", "N" = "N", "P" = "P", "Q" = "Q", "R" = "R", "S" = "S", "T" = "T", "V" = "V", "W" = "W", "Y" = "Y"))
    ''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "MutantAA", title = title, plot_scale = plot_scale, point_opacity = 0.6)


    def scatterplot_GP(self, title, csv_filename):

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'HasGPMutation']]
        new_dataframe['GP'] = numpy.where(new_dataframe['HasGPMutation'] == 1, 'GP', 'Other')
        new_dataframe['Opacity'] = numpy.where(new_dataframe['HasGPMutation'] == 1, 0.9, 0.5)
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'GP', 'Opacity']]
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''plot_scale <- scale_color_manual(
    name="Glycine/Proline",
    values = c( "None" = '#777777', "GP" = '%(neon_green)s', "Other" = '#440077'),
    labels = c( "None" = "N/A", "GP" = "GP", "Other" = "Other"))''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "GP", title = title, plot_scale = plot_scale, point_opacity = 0.75)


    def scatterplot_pdb_res_binned(self, title, csv_filename):
        '''Scatterplot by binned PDB resolution.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'PDBResolutionBin']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    name = "Resolution",
    values = c( "N/A" = '#777777', "<1.5" = '#0052aE', "1.5-2.0" = '#554C54', '2.0-2.5' = "#FFA17F", '>=2.5' = "#ce4200")
    )''' % plot_colors
        return self.scatterplot_color_by_series(colorseries = "PDBResolutionBin", title = title, plot_scale = plot_scale, point_opacity = 0.75)


    def scatterplot_chain_length(self, title, csv_filename):
        '''Scatterplot by chain length.'''

        # Create CSV input
        new_dataframe = self.dataframe.copy()
        new_dataframe = new_dataframe[['Experimental', 'Predicted', 'MonomerLength']]
        new_dataframe['Opacity'] = 0.4
        new_dataframe.columns = ['Experimental', 'Predicted', 'Residues', 'Opacity'] # rename the monomer length column
        new_dataframe.to_csv(csv_filename, sep = ',', header = True)

        plot_scale = '''
plot_scale <- scale_color_manual(
    name = "Resolution",
    values = c( "N/A" = '#777777', "<1.5" = '#0052aE', "1.5-2.0" = '#554C54', '2.0-2.5' = "#FFA17F", '>=2.5' = "#ce4200")
    )''' % plot_colors
        extra_commands ='\n    scale_colour_gradient(low="yellow", high="#880000") +'
        return self.scatterplot_color_by_series(colorseries = "Residues", title = title, plot_scale = '', point_opacity = 0.75, extra_commands = extra_commands)
