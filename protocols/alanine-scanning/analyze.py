#!/usr/bin/python

import sys
import os
import argparse
import cPickle as pickle
import sqlite3
import itertools

import analysis.stats as stats
from alanine_scanning import parse_mutations_file

pickle_name = os.path.join('data', 'job_dict.pickle')
global_analysis_output_dir = 'analysis_output'

def parse_db_output(db_output_file, ddg_data, score_fxns):
    query = 'SELECT ddG.resNum, ddg.mutated_to_name3, ddg.ddG_value, ' \
    'structures.tag, residue_pdb_identification.residue_number, ' \
    'residue_pdb_identification.chain_id, residues.name3, batches.description, batches.name ' \
    'FROM ddg INNER JOIN structures ON structures.struct_id=ddg.struct_id ' \
    'INNER JOIN residue_pdb_identification ON residue_pdb_identification.struct_id=structures.struct_id ' \
    'AND residue_pdb_identification.residue_number=ddg.resNum ' \
    'INNER JOIN residues ON residues.struct_id=structures.struct_id AND residues.resNum=ddg.resNum ' \
    'INNER JOIN batches ON batches.batch_id=structures.batch_id'

    conn = sqlite3.connect(db_output_file)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    for rosetta_resnum, mutated_to, ddg_calc, tag, pdb_resnum, chain, original_name3, description, PDBPosID in c.execute(query):
        pdb_id = tag.split('_')[0]
        assert( len(pdb_id) == 4 )
        if PDBPosID not in ddg_data:
            ddg_data[PDBPosID] = {}
        ddg_data[PDBPosID][description] = ddg_calc
        score_fxns.add(description)

    conn.close()
    return (ddg_data, score_fxns)

def get_mutations_data_from_pdb_dict(pdb_data_dict, PDBPosID):
    for m in pdb_data_dict.values():
        if PDBPosID in m.PDBPosID_list:
            return m
    raise Exception("Couldn't match ID %s" % str(PDBPosID))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('output_dirs',
                        nargs='*',
                        help = 'Output directories')
    args = parser.parse_args()

    mutations_data = parse_mutations_file()

    for output_dir in args.output_dirs:

        analysis_output_dir = os.path.join(global_analysis_output_dir, os.path.basename(output_dir))
        if not os.path.isdir(analysis_output_dir):
            os.makedirs(analysis_output_dir)

        assert( os.path.isdir(output_dir) )
        job_pickle_file = os.path.join(
            output_dir, pickle_name
        )
        
        with open(job_pickle_file,'r') as p:
            job_dict = pickle.load(p)

        db_output_files = []
        for job_dir in job_dict:
            job_dir = os.path.join(output_dir, job_dir)
            db_output_file = os.path.join(job_dir, 'rosetta_output.db3')
            if os.path.isfile(db_output_file):
                db_output_files.append(db_output_file)

        print 'Found %d output dbs (%d expected)' % (len(db_output_files), len(job_dict))

        score_fxns = set()
        ddg_data = {}
        for db_output_file in db_output_files:
            ddg_data, score_fxns = parse_db_output(db_output_file, ddg_data, score_fxns)

        for mut in mutations_data.values():
            for i, PDBPosID in enumerate(mut.PDBPosID_list):
                if PDBPosID not in ddg_data:
                    ddg_data[PDBPosID] = {}
                ddg_data[PDBPosID]['ddg_obs'] = mut.ddg_obs_list[i]
                ddg_data[PDBPosID]['tanja_ddg_calc'] = mut.ddg_calc_list[i]

        score_fxns = sorted(list(score_fxns))
        score_fxns.insert(0, 'tanja_ddg_calc')
        score_fxns.insert(0, 'ddg_obs')

        # Write output CSV and create data lists for further analysis
        all_data_ids = []
        data_id_in_interface = [] # True/False mask array
        all_data_points = [[] for x in xrange(len(score_fxns))]
        with open(
            os.path.join(analysis_output_dir,
                         'results-%s.csv' % os.path.basename(output_dir)),
            'w'
        ) as f:
            f.write('ID')
            for score_fxn in score_fxns:
                f.write(',%s' % score_fxn)
            f.write('\n')
            for PDBPosID in sorted( ddg_data.keys() ):
                f.write('%s' % PDBPosID)
                data_is_complete = True
                for score_fxn in score_fxns:
                    if score_fxn in ddg_data[PDBPosID]:
                        f.write(',%.6f' % ddg_data[PDBPosID][score_fxn])
                    else:
                        data_is_complete = False
                        f.write(',NA')
                if data_is_complete:
                    all_data_ids.append(PDBPosID)
                    for i, score_fxn in enumerate(score_fxns):
                        all_data_points[i].append(ddg_data[PDBPosID][score_fxn])
                    if get_mutations_data_from_pdb_dict(mutations_data, PDBPosID).in_interface_by_id(PDBPosID):
                        data_id_in_interface.append(True)
                    else:
                        data_id_in_interface.append(False)

                f.write('\n')

        # Run stats module helper
        def run_stats(name, i_name, j_name, data_ids, i_data_points, j_data_points):
            dataset_stats = stats._get_xy_dataset_statistics(
                i_data_points, j_data_points
            )

            stats_str = stats.format_stats_for_printing(dataset_stats)
            with open(os.path.join(analysis_output_dir, '%s-stats.txt' % os.path.basename(output_dir)), 'a') as f:
                f.write('%s - %s vs %s\n' % (name, i_score_fxn, j_score_fxn) )
                f.write(stats_str)
                f.write('\n\n')
            print stats_str
            print

            table_for_plot = []
            for pt_id, pt_i, pt_j in zip(data_ids, i_data_points, j_data_points):
                table_for_plot.append(dict(ID = pt_id, Experimental = pt_i, Predicted = pt_j))

            stats.plot(
                table_for_plot,
                os.path.join(
                    analysis_output_dir,
                    '%s-%s-%s_vs_%s.pdf' % ( os.path.basename(output_dir), name, i_name, j_name)
                ),
                stats.RInterface.correlation_coefficient_gplot
            )

        def compress(data, selectors, invert=False):
            # Compress implementation with inversion option
            # compress('ABCDEF', [1,0,1,0,1,1]) --> A C E F
            if invert:
                return [d for d, s in itertools.izip(data, selectors) if not s]
            else:
                return [d for d, s in itertools.izip(data, selectors) if s]

        # Uncomment this outer loop to do all-by-all comparison
        # Otherwise, default is to only compare against first, experimental results
        # for i, i_score_fxn in enumerate(score_fxns):
        i = 0
        i_score_fxn = score_fxns[0]
        for j, j_score_fxn in enumerate(score_fxns):
            if i == j:
                continue

            print '########', i_score_fxn, 'vs', j_score_fxn, '########'

            print '#### Points in interface: ####'
            run_stats(
                'interface_pts', score_fxns[i], score_fxns[j],
                compress(all_data_ids, data_id_in_interface),
                compress(all_data_points[i], data_id_in_interface),
                compress(all_data_points[j], data_id_in_interface)
            )

            # print '#### Points not in interface: ####'
            # run_stats(
            #     'noninterface_pts', score_fxns[i], score_fxns[j],
            #     compress(all_data_ids, data_id_in_interface, invert=True),
            #     compress(all_data_points[i], data_id_in_interface, invert=True),
            #     compress(all_data_points[j], data_id_in_interface, invert=True)
            # )

            # print '#### All points: ####'
            # run_stats(
            #     'all_pts', score_fxns[i], score_fxns[j],
            #     all_data_ids,
            #     all_data_points[i],
            #     all_data_points[j]
            # )
