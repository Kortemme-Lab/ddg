#!/usr/bin/python

import multiprocessing
import os
from rosetta.write_run_file import process as write_run_file
import shutil
import sys
import time
import inspect
import cPickle as pickle
import re
import getpass
import interfaces_defs

input_pdb_dir_path = '../../input/pdbs/hydrogen_pdbs'
extra_name = 'rp_rp-bound_rlx-bound-cart' # something like _talaris if needed
mutations_file_location = 'MUTATIONS.dat'
rosetta_scripts_protocol = 'alascan.xml'
resfile_start = 'NATRO\nEX 1 EX 2 EX 3\nSTART\n'
score_fxns = ['talaris2014', 'soft_rep', 'talaris2013', 'talaris2014_soft_fa_rep', 'talaris2014_soft_fa_dun', 'score12', 'interface']
cluster_rosetta_bin = '/netapp/home/kbarlow/alascan/source/bin'
local_rosetta_bin = '/home/kyleb/rosetta/working_branches/alascan/source/bin'
job_output_directory = 'job_output'

class MutationData:
    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        self.chainres_list = []
        self.pdb_res_list = []
        self.amino_acid_list = []
        self.ddg_calc_list = []
        self.ddg_obs_list = []
        self.ddg_obs_greater_than_list = []
        self.tanja_id_list = []
        self.chain_list = []
        self.interface_list = []
        self.new_id_list = []
        self.insertion_code_list = []
        self.id_conv = {}

        for pdb_dict in interfaces_defs.kortemme_baker_protein_protein_interfaces:
            if pdb_dict['PDBFileID'] == pdb_id:
                self.lchains = set(pdb_dict['LChains'])
                self.rchains = set(pdb_dict['RChains'])
                break

        # Load chain order
        self.chain_order = get_chain_order(pdb_id)

    def tanja_id_to_alascan(self, tanja_id):
        i = self.tanja_id_list.index(tanja_id)
        # print self.chainres_list[i][1:]
        # print self.pdb_res_list[i]
        return ( three_letter_codes[self.amino_acid_list[i]], int(self.pdb_res_list[i]), self.chain_list[i] )

    def add_line(self, line):
        data = line.split()
        pdb_id = data[0]
        assert( pdb_id == self.pdb_id )
        chain = data[2]

        if data[3] == '100B':
            pdb_res = 100
            insertion_code = 'B'
        else:
            pdb_res = int(data[3])
            insertion_code = ''

        amino_acid = data[4]
        assert( len(amino_acid) == 1 )
        ddg_calc = float(data[5])

        if data[6].startswith('>'):
            greater_than = True
            ddg_obs = float( data[6][1:] )
        else:
            greater_than = False
            ddg_obs = float(data[6])

        if data[7] == '0':
            interface = False
        elif data[7] == '1':
            interface = True

        tanja_id = data[8]
        assert( tanja_id.startswith(pdb_id) )

        new_id = '%s%d' % (pdb_id, pdb_res)

        chainres = '%s%d' % (chain, pdb_res)
        assert(chainres not in self.chainres_list)
        self.chainres_list.append(chainres)

        self.pdb_res_list.append(pdb_res)
        self.amino_acid_list.append(amino_acid)
        self.ddg_calc_list.append(ddg_calc)
        self.ddg_obs_list.append(ddg_obs)
        self.ddg_obs_greater_than_list.append(greater_than)
        self.tanja_id_list.append(tanja_id)
        self.chain_list.append(chain)
        self.interface_list.append(interface)
        self.new_id_list.append(new_id)
        self.insertion_code_list.append(insertion_code)
        self.id_conv[tanja_id] = new_id

    def in_interface(self, resnum):
        return self.interface_list[self.pdb_res_list.index(resnum)]

    def num_chains(self):
        return len( set(self.chain_list) )

    def num_residues(self):
        return len(self.pdb_res_list)

    def lchainsnumbers(self):
        # Figure out lchains order
        lchainsnumbers = ''
        for lchain in self.lchains:
            lchainsnumbers += str(self.chain_order.index(lchain)+1)
        return lchainsnumbers

    def rchainsnumbers(self):
        # Figure out rchains order
        rchainsnumbers = ''
        for rchain in self.rchains:
            rchainsnumbers += str(self.chain_order.index(rchain)+1)
        return rchainsnumbers

    def get_comma_chains_to_move(self):
        lchainsnumbers = self.lchainsnumbers()
        rchainsnumbers = self.rchainsnumbers()

        if '1' in lchainsnumbers:
            assert( '1' not in rchainsnumbers )
            chains_to_move = rchainsnumbers
        else:
            chains_to_move = lchainsnumbers

        comma_chains_to_move = ''
        for i, char in enumerate(chains_to_move):
            comma_chains_to_move += char
            if i+1 < len(chains_to_move):
                comma_chains_to_move += ','
        
        return comma_chains_to_move

    def __repr__(self):
        return '%s: %d chains, %d residues' % (self.pdb_id, self.num_chains(), self.num_residues())

def get_chain_order(pdb_id):
    allchains_pdb = os.path.join(input_pdb_dir_path, '%s_0001.pdb' % pdb_id)
    last_chain = None
    chain_order = ''
    all_chains_set = set()
    with open(allchains_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                all_chains_set.add(chain)
                if chain != last_chain:
                    last_chain = chain
                    chain_order += chain

    assert( len(all_chains_set) == len(chain_order) )
    return chain_order

def parse_mutations_file():
    return_dict = {}
    with open(mutations_file_location, 'r') as f:
        for line in f:
            pdb_id = line.split()[0]
            if pdb_id == 'PDB_ID': # Skip first line
                continue
            if pdb_id not in return_dict:
                return_dict[pdb_id] = MutationData(pdb_id)
            return_dict[pdb_id].add_line(line)
    return return_dict

def make_resfile(resfile_path, mutation_datum, tanja_id):
    # Make resfile
    mutation_info_index = mutation_datum.tanja_id_list.index(tanja_id)
    chain = mutation_datum.chain_list[mutation_info_index]
    pdb_res = mutation_datum.pdb_res_list[mutation_info_index]
    insertion_code = mutation_datum.insertion_code_list[mutation_info_index]

    with open(resfile_path, 'w') as f:
        f.write(resfile_start)
        f.write('%d%s %s PIKAA A\n' % (pdb_res, insertion_code, chain))

if __name__ == "__main__":
    mutation_info = parse_mutations_file()

    job_name = os.path.basename(inspect.getfile(inspect.currentframe())).split('.')[0] + extra_name
    output_dir = os.path.join(job_output_directory, '%s-%s_%s' % (time.strftime("%y%m%d"), getpass.getuser(), job_name) )
    output_data_dir = os.path.join(output_dir, 'data')
    pdb_data_dir = os.path.join(output_data_dir, 'input_pdbs')
    
    if not os.path.isdir(pdb_data_dir):
        os.makedirs(pdb_data_dir)

    # Copy Rosetta scripts protocol
    protocol_path = os.path.join(output_data_dir, os.path.basename(rosetta_scripts_protocol))
    shutil.copy(rosetta_scripts_protocol, protocol_path)
    protocol_relpath = os.path.relpath(protocol_path, output_dir)

    resfile_data_dir = os.path.join(output_data_dir, 'resfiles')
    if not os.path.isdir(resfile_data_dir):
        os.makedirs(resfile_data_dir)

    job_dict = {}
    for pdb_id in mutation_info:
        pdb_path = os.path.join(input_pdb_dir_path, '%s_0001.pdb' % pdb_id)
        new_pdb_path = os.path.join(pdb_data_dir, os.path.basename(pdb_path))
        if not os.path.isfile(new_pdb_path):
            # Assume file is correct version to save copy time
            shutil.copy(pdb_path, new_pdb_path)
        pdb_relpath = os.path.relpath(new_pdb_path, output_dir)

        # Setup chains to move
        comma_chains_to_move = mutation_info[pdb_id].get_comma_chains_to_move()
                
        for score_fxn in score_fxns:    
            for tanja_id in mutation_info[pdb_id].tanja_id_list:
                # Make resfile
                resfile_path = os.path.join(resfile_data_dir, '%s.mutation.resfile' % tanja_id)
                if not os.path.exists(resfile_path):
                    make_resfile(resfile_path, mutation_info[pdb_id], tanja_id)
                resfile_relpath = os.path.relpath(resfile_path, output_dir)

                sub_dict = {}
                sub_dict['-in:file:s'] = pdb_relpath

                sub_dict['-parser:protocol'] = protocol_relpath

                sub_dict['-parser:script_vars'] = [
                    'tanja_id=%s' % tanja_id,
                    'currentscorefxn=%s' % score_fxn,
                    'chainstomove=%s' % comma_chains_to_move,
                    'currentpackscorefxn=%s' % score_fxn,
                    'pathtoresfile=%s' % resfile_relpath,
                ]

                job_dict[ '%s/%s/%s' % (pdb_id.upper(), tanja_id, score_fxn) ] = sub_dict

    with open(os.path.join(output_data_dir, 'job_dict.pickle'), 'w') as f:
        pickle.dump(job_dict, f)

    general_rosetta_args = "'-parser:view', '-inout:dbms:mode', 'sqlite3', '-inout:dbms:database_name', 'rosetta_output.db3', '-no_optH', 'true'"

    args = {
        'scriptname' : 'ddg_run',
        'appname' : 'rosetta_scripts.mysql.linuxgccrelease',
        'rosetta_args_list' : general_rosetta_args,
    }

    args['cluster_rosetta_bin'] = cluster_rosetta_bin
    args['local_rosetta_bin'] = local_rosetta_bin
    args['numjobs'] = '%d' % len(job_dict)
    args['mem_free'] = '1.1G'
    args['output_dir'] = output_dir

    write_run_file(args)

    print 'Job files written to directory:', os.path.abspath(output_dir)
