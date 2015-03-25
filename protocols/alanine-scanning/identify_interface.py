#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Kyle A. Barlow
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

import os
import sys
import setup_alanine_scanning
import cPickle as pickle
from math import sqrt

cutoff_dist = 5.0

pdb_atom_dict = {}

def coords_from_line(line):
    '''
    Returns a tuple (x,y,z) of floats of coordinates from a string
    line from a pdb
    '''
    return ((
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
        ))

def distance_3d(coord1,coord2):
    '''
    Returns the euclidean distance between (x,y,z) tuple pairs coord1 and coord2
    '''
    return sqrt( (coord1[0]-coord2[0])**2 +
                 (coord1[1]-coord2[1])**2 +
                 (coord1[2]-coord2[2])**2 )

def read_pdb_atoms(pdb_file):
    atoms_list = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                element = line[76:78].strip()
                if element == 'H':
                    continue
                chain = line[21]
                resnum = int(line[22:26].strip())
                insertion_code = line[26].strip()
                coords = coords_from_line(line)
                atoms_list.append( (chain, resnum, insertion_code, coords) )

    return atoms_list

def find_resnums_around_mutation(pdb_id, resnum, insertion_code, chain):
    if pdb_id not in pdb_atom_dict:
        pdb_atom_dict[pdb_id] = read_pdb_atoms( os.path.join(setup_alanine_scanning.input_pdb_dir_path, '%s_0001.pdb' % pdb_id) )
    pdb_atoms = pdb_atom_dict[pdb_id]

    mut_coords = []
    for atom_chain, atom_resnum, atom_insertion_code, atom_coords in pdb_atoms:
        if atom_resnum == resnum and insertion_code == atom_insertion_code and chain == atom_chain:
            mut_coords.append( ((atom_resnum, atom_insertion_code, atom_chain), atom_coords) )

    assert( len(mut_coords) != 0 )

    def check_close_helper(atom_chain, atom_resnum, atom_insertion_code, atom_coords):
        for info_tuple, mut_coord in mut_coords:
            resnum, insertion_code, chain = info_tuple
            if atom_resnum == resnum and insertion_code == atom_insertion_code and chain == atom_chain:
                continue
            if distance_3d(mut_coord, atom_coords) <= cutoff_dist:
                return (atom_resnum, atom_insertion_code, atom_chain)
        return None

    close_atoms = []
    close_resnums = set()
    for atom_chain, atom_resnum, atom_insertion_code, atom_coords in pdb_atoms:
        if atom_resnum in close_resnums:
            continue
        result = check_close_helper(atom_chain, atom_resnum, atom_insertion_code, atom_coords)
        if result != None:
            close_atoms.append(result)
            close_resnums.add(atom_resnum)

    return list(close_atoms)

def get_close_residues_dict():
    if os.path.isfile('.close_residues.pickle'):
        with open('.close_residues.pickle', 'r') as f:
            return pickle.load(f)

    mut_data_dict = setup_alanine_scanning.parse_mutations_file()
    close_residues_dict = {}
    for pdb_id in mut_data_dict:
        mut_data = mut_data_dict[pdb_id]

        for i, pdb_res in enumerate(mut_data.pdb_res_list):
            insertion_code = mut_data.insertion_code_list[i]
            chain = mut_data.chain_list[i]
            close_residues = find_resnums_around_mutation(pdb_id, int(pdb_res), insertion_code, chain)

            print 'Mutation %s:%d%s has %d residues within %.1f angstroms' % (pdb_id, pdb_res, insertion_code, len(close_residues), cutoff_dist)

            if pdb_id not in close_residues_dict:
                close_residues_dict[pdb_id] = {}
            reskey = (pdb_res, chain, insertion_code)
            assert( reskey not in close_residues_dict[pdb_id] )
            close_residues_dict[pdb_id][reskey] = close_residues

    with open('.close_residues.pickle', 'w') as f:
        pickle.dump(close_residues_dict, f)

    return close_residues_dict
