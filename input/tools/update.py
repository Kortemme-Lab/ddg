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


import os
import json

from klab.bio.pdb import PDB
from klab.fs.fsio import read_file, write_file


def update_pdbs_json():
    '''This function was used to update the pdbs.json file to include chain sequences and types.'''
    pdb_data = {}
    pdb_data_ = json.loads(read_file(os.path.join('..', 'json', 'pdbs.json')))
    for k, v in pdb_data_.iteritems():
        assert(len(k) == 4)
        newk = k.upper()
        pdb = PDB(read_file(os.path.join('..', 'pdbs', newk + '.pdb')))
        chain_ids = set(pdb.chain_types.keys()).union(set(pdb.seqres_chain_order)).union(set(pdb.atom_sequences.keys()))
        v['Chains'] = dict.fromkeys(chain_ids)
        for chain_id in chain_ids:
            v['Chains'][chain_id] = dict(
                Sequence = str(pdb.atom_sequences.get(chain_id)),
                Type = pdb.chain_types.get(chain_id),
            )
        pdb_data[newk] = v
    write_file(os.path.join('..', 'json', 'pdbs.json.new'), json.dumps(pdb_data, indent = 4, sort_keys=True))

