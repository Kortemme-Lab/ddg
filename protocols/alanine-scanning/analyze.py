#!/usr/bin/python

import sys
import os
import argparse
import cPickle as pickle
import sqlite3

pickle_name = os.path.join('data', 'job_dict.pickle')

def parse_db_output(db_output_file, ddg_data):
    query = 'SELECT ddG.resNum, ddg.mutated_to_name3, ddg.ddG_value, ' \
    'structures.tag, residue_pdb_identification.residue_number, ' \
    'residue_pdb_identification.chain_id, residues.name3, batches.description ' \
    'FROM ddg INNER JOIN structures ON structures.struct_id=ddg.struct_id ' \
    'INNER JOIN residue_pdb_identification ON residue_pdb_identification.struct_id=structures.struct_id ' \
    'AND residue_pdb_identification.residue_number=ddg.resNum ' \
    'INNER JOIN residues ON residues.struct_id=structures.struct_id AND residues.resNum=ddg.resNum ' \
    'INNER JOIN batches ON batches.batch_id=structures.batch_id'

    conn = sqlite3.connect(db_output_file)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    for rosetta_resnum, mutated_to, ddg_calc, tag, pdb_resnum, chain, original_name3, description in c.execute(query):
        pdb_id = tag.split('_')[0]
        assert( len(pdb_id) == 4 )
        ddg_data[(pdb_id, original_name3, pdb_resnum, chain, mutated_to, pdb_resnum, description)] = ddg_calc

    conn.close()
    return ddg_data

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('output_dirs',
                        nargs='*',
                        help = 'Output directories')
    args = parser.parse_args()

    ddg_data = {}
    for output_dir in args.output_dirs:
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

        for db_output_file in db_output_files:
            ddg_data = parse_db_output(db_output_file, ddg_data)

    print 'Processed %d data points' % len(ddg_data)
    with open('results.csv', 'w') as f:
        f.write('PDB,Mutation,score_function,ddg_calc\n')
        for tup in ddg_data:
            pdb_id, original_name3, pdb_resnum, chain, mutated_to, pdb_resnum, description = tup
            f.write('%s%d,%s%d%s%s,%s,%.6f\n' % (pdb_id, pdb_resnum, original_name3, pdb_resnum, chain, mutated_to, description, ddg_data[tup]) )
