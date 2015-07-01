#!/usr/bin/env python2
# encoding: utf-8

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

import re
import sys
import os
import types
import string
import types

from basics import Residue, PDBResidue, Sequence, SequenceMap, residue_type_3to1_map, protonated_residue_type_3to1_map, non_canonical_amino_acids, protonated_residues_types_3, residue_types_3, Mutation, ChainMutation, SimpleMutation
from basics import dna_nucleotides, rna_nucleotides, dna_nucleotides_3to1_map, dna_nucleotides_2to1_map, non_canonical_dna, non_canonical_rna, all_recognized_dna, all_recognized_rna
from analysis.stats import read_file, write_file
from map_pdb_residues import get_pdb_contents_to_pose_residue_map


### Residue types

allowed_PDB_residues_types = protonated_residues_types_3.union(residue_types_3)
allowed_PDB_residues_and_nucleotides = allowed_PDB_residues_types.union(dna_nucleotides).union(rna_nucleotides)

### Rosetta patches

# Some Rosetta versions fails on some edge cases with certain residues. Since we rely on a lot of logic associated with this module (mappings between residues), it seems best to fix those here.
ROSETTA_HACKS_residues_to_remove = {
    '1A2P' : set(['B   3 ']), # terminal residue B 3 gets removed which triggers an exception "( anchor_rsd.is_polymer() && !anchor_rsd.is_upper_terminus() ) && ( new_rsd.is_polymer() && !new_rsd.is_lower_terminus() ), Exit from: src/core/conformation/Conformation.cc line: 449". This seems like a Rosetta deficiency.
    '1BAO' : set(['B   3 ']), # similar
    '1BNI' : set(['C   3 ']), # similar
    '1BRH' : set(['A   3 ', 'B   3 ', 'C   3 ']), # similar
    '1BRI' : set(['A   3 ', 'B   3 ']), # similar
    '1BRJ' : set(['A   3 ', 'B   3 ', 'C   3 ']), # similar
    '1BRK' : set(['B   3 ']), # similar
    '1EHK' : set(['B   3 ']), # similar
    '1ZNJ' : set(['B  30 ', 'F   1 ']), # similar
    '487D' : set(['N   1 ']),

}

### PDB hacks

missing_chain_ids = {
    '2MBP' : 'A', # The FASTA file lists this as 'A' so we need to patch records up to match
}

### Whitelist for PDB files with ACE residues (we could allow all to pass but it may be good to manually look at each case)

cases_with_ACE_residues_we_can_ignore = set(['3UB5', '1TIN', '2ZTA', '5CPV', '1ATN', '1LFO', '1OVA', '3PGK', '2FAL', '2SOD', '1SPD', '1UOX', '1UNQ', '1DFJ', '1B39', '1HRC', '1PNE', '1WEJ', '2BNH'])


### Some PDB files have been deprecated but not replaced. Again, it may be good to manually look at each case.

obsolete_pdb_ids_with_no_replacement_entries = set(['1CMW'])


### Parsing-related variables

COMPND_field_map = {
    'MOL_ID' : 'MoleculeID',
    'MOLECULE' : 'Name',
    'CHAIN' : 'Chains',
    'FRAGMENT' : 'Fragment',
    'SYNONYM' : 'Synonym',
    'EC' : 'EC',
    'ENGINEERED' : 'Engineered',
    'MUTATION' : 'Mutation',
    'OTHER_DETAILS' : 'OtherDetails',
}

SOURCE_field_map = {
    'MOL_ID' : 'MoleculeID',
    'SYNTHETIC' : 'Synthetic',
    'ORGANISM_SCIENTIFIC' : 'OrganismScientificName',
    'ORGANISM_COMMON' : 'OrganismCommonName',
    'ORGANISM_TAXID' : 'OrganismNCBITaxonomyID',
}

modified_residues_patch = {
    '1A2C' : {
        '34H' : 'UNK',
    },
    '2ATC' : {
        'ASX' : 'ASN',
    },
    '1XY1' : {
        'MPT' : 'UNK',
    },
    '1CVW' : { # Note: more recent versions of this file do not require this patch
        'ANS' : 'UNK',
        '0QE' : 'UNK',
    },
    '1FAK' : {
        'CGU' : 'GLU', # Gamma-carboxy-glutamic acid
    },
    '1JXQ' : {
        'PHQ' : 'UNK', # benzyl chlorocarbonate
        'CF0' : 'UNK', # fluoromethane
    },
    '1YJ1' : {
        'DGN' : 'GLN', # D-glutamine
    },
    '2CN0' : {
        'SIN' : 'UNK', # Succinic acid
    },
    '2FTL' : {
        'IAS' : 'ASP', # Beta-L-aspartic acid/L-isoaspartate. Mismatch to asparagine - "the expected l-Asn residue had been replaced with a non-standard amino acid" (10.1016/j.jmb.2006.11.003).
    },
}

### Record types

order_of_records = [
    "HEADER","OBSLTE","TITLE","SPLIT","CAVEAT","COMPND","SOURCE","KEYWDS",
    "EXPDTA","NUMMDL","MDLTYP","AUTHOR","REVDAT","SPRSDE","JRNL","REMARK",
    "DBREF","DBREF1","DBREF2","DBREF1/DBREF2","SEQADV","SEQRES","MODRES",
    "HET","HETNAM","HETSYN","FORMUL","HELIX","SHEET","SSBOND","LINK","CISPEP",
    "SITE","CRYST1","ORIGX1","ORIGX2","ORIGX3","SCALE1","SCALE2","SCALE3",
    "MTRIX1","MTRIX2","MTRIX3","MODEL","ATOM","ANISOU","TER","HETATM",
    "ENDMDL","CONECT","MASTER","END"
]
order_of_records = [x.ljust(6) for x in order_of_records]

allowed_record_types = set([
# One time, single line:
'CRYST1', #     Unit cell parameters, space group, and Z.
'END   ', #     Last record in the file.
'HEADER', #     First line of the entry, contains PDB ID code, classification, and date of deposition.
'NUMMDL', #     Number of models.
'MASTER', #     Control record for bookkeeping.
'ORIGXn', #     Transformation from orthogonal  coordinates to the submitted coordinates (n = 1, 2, or 3).
'SCALEn', #     Transformation from orthogonal coordinates to fractional crystallographic coordinates  (n = 1, 2, or 3).
# One time, multiple lines:
'AUTHOR', #     List of contributors.
'CAVEAT', #     Severe error indicator.
'COMPND', #     Description of macromolecular contents of the entry.
'EXPDTA', #     Experimental technique used for the structure determination.
'MDLTYP', #     Contains additional annotation  pertinent to the coordinates presented  in the entry.
'KEYWDS', #     List of keywords describing the macromolecule.
'OBSLTE', #     Statement that the entry has been removed from distribution and list of the ID code(s) which replaced it.
'SOURCE', #     Biological source of macromolecules in the entry.
'SPLIT ', #     List of PDB entries that compose a larger  macromolecular complexes.
'SPRSDE', #     List of entries obsoleted from public release and replaced by current entry.
'TITLE ', #     Description of the experiment represented in the entry.
# Multiple times, one line:
'ANISOU', #     Anisotropic temperature factors.
'ATOM  ', #     Atomic coordinate records for  standard groups.
'CISPEP', #     Identification of peptide residues in cis conformation.
'CONECT', #     Connectivity records.
'DBREF ', #     Reference  to the entry in the sequence database(s).
'HELIX ', #     Identification of helical substructures.
'HET   ', #     Identification of non-standard groups heterogens).
'HETATM', #     Atomic coordinate records for heterogens.
'LINK  ', #     Identification of inter-residue bonds.
'MODRES', #     Identification of modifications to standard residues.
'MTRIXn', #     Transformations expressing non-crystallographic symmetry (n = 1, 2, or 3). There may be multiple sets of these records.
'REVDAT', #     Revision date and related information.
'SEQADV', #     Identification of conflicts between PDB and the named sequence database.
'SHEET ', #     Identification of sheet substructures.
'SSBOND', #     Identification of disulfide bonds.
# Multiple times, multiple lines:
'FORMUL', #     Chemical formula of non-standard groups.
'HETNAM', #     Compound name of the heterogens.
'HETSYN', #     Synonymous compound names for heterogens.
'SEQRES', #     Primary sequence of backbone residues.
'SITE  ', #     Identification of groups comprising important entity sites.
# Grouping:
'ENDMDL', #     End-of-model record for multiple structures in a single coordinate entry.
'MODEL ', #     Specification of model number for multiple structures in a single coordinate entry.
'TER   ', #     Chain terminator.
# Other:
'JRNL  ', #     Literature citation that defines the coordinate set.
'REMARK', #     General remarks; they can be structured or free form.
])

# This set is probably safer to use to allow backwards compatibility
all_record_types = allowed_record_types.union(set(order_of_records))

### Exception classes
class PDBParsingException(Exception): pass
class MissingRecordsException(Exception): pass
class NonCanonicalResidueException(Exception): pass
class PDBValidationException(Exception): pass
class PDBMissingMainchainAtomsException(Exception): pass


def create_mutfile(pdb, mutations):
    '''The mutations here are in the original PDB numbering. The PDB file will use Rosetta numbering.
        We use the mapping from PDB numbering to Rosetta numbering to generate the mutfile.
    '''
    mutfile = []

    for mutation in mutations:
        chain = mutation.Chain
        resid = PDB.ResidueID2String(mutation.ResidueID)
        wt = mutation.WildTypeAA
        mt = mutation.MutantAA

        # Check that the expected wildtype exists in the PDB
        readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
        assert(wt == residue_type_3to1_map[readwt])
        resid = resid.strip()
        mutfile.append("%(wt)s %(resid)s %(mt)s" % vars())
    if mutfile:
        mutfile = ["total %d" % len(mutations), "%d" % len(mutations)] + mutfile
        return '\n'.join(mutfile)
    else:
        raise Exception("An error occurred creating the mutfile.")


class PDB:
    """A class to store and manipulate PDB data"""

    ### Constructor ###

    def __init__(self, pdb_content, pdb_id = None, strict = True):
        '''Takes either a pdb file, a list of strings = lines of a pdb file, or another object.'''

        self.pdb_content = pdb_content
        if type(pdb_content) is types.StringType:
            self.lines =  pdb_content.split("\n")
        else:
            self.lines = [line.strip() for line in pdb_content]
        self.parsed_lines = {}
        self.structure_lines = []                       # For ATOM and HETATM records
        self.journal = None
        self.chain_types = {}
        self.format_version = None
        self.modified_residues = None
        self.modified_residue_mapping_3 = {}
        self.pdb_id = None
        self.strict = strict

        self.seqres_chain_order = []                    # A list of the PDB chains in document-order of SEQRES records
        self.seqres_sequences = {}                      # A dict mapping chain IDs to SEQRES Sequence objects
        self.atom_chain_order = []                      # A list of the PDB chains in document-order of ATOM records (not necessarily the same as seqres_chain_order)
        self.atom_sequences = {}                        # A dict mapping chain IDs to ATOM Sequence objects
        self.chain_atoms = {}                           # A dict mapping chain IDs to a set of ATOM types. This is useful to test whether some chains only have CA entries e.g. in 1LRP, 1AIN, 1C53, 1HIO, 1XAS, 2TMA

        # PDB deprecation fields
        self.deprecation_date = None
        self.deprecated = False
        self.replacement_pdb_id = None

        self.rosetta_to_atom_sequence_maps = {}
        self.rosetta_residues = []
        self.residue_types = set()                      # the set of 3-letter residue types present in the file (useful for checking against e.g. CSE, MSE)

        self.fix_pdb()
        self._split_lines()
        self.pdb_id = pdb_id
        self.pdb_id = self.get_pdb_id()                 # parse the PDB ID if it is not passed in
        self._apply_hacks()
        self._get_pdb_format_version()
        self._get_modified_residues()
        self._get_replacement_pdb_id()
        if missing_chain_ids.get(self.pdb_id):
            self._update_structure_lines()
        self._get_SEQRES_sequences()
        self._get_ATOM_sequences()


    def fix_pdb(self):
        '''A function to fix fatal errors in PDB files when they can be automatically fixed. At present, this only runs if
           self.strict is False. We may want a separate property for this since we may want to keep strict mode but still
           allow PDBs to be fixed.

           The only fixes at the moment are for missing chain IDs which get filled in with a valid PDB ID, if possible.'''

        if self.strict:
            return

        # Get the list of chains
        chains = set()
        for l in self.lines:
            if l.startswith('ATOM  ') or l.startswith('HETATM'):
                chains.add(l[21])

        # If there is a chain with a blank ID, change that ID to a valid unused ID
        if ' ' in chains:
            fresh_id = None
            allowed_chain_ids = list(string.uppercase) + list(string.lowercase) + map(str, range(10))
            for c in chains:
                try: allowed_chain_ids.remove(c)
                except: pass
            if allowed_chain_ids:
                fresh_id = allowed_chain_ids[0]

            # Rewrite the lines
            new_lines = []
            if fresh_id:
                for l in self.lines:
                    if (l.startswith('ATOM  ') or l.startswith('HETATM')) and l[21] == ' ':
                        new_lines.append('%s%s%s' % (l[:21], fresh_id, l[22:]))
                    else:
                        new_lines.append(l)
                self.lines = new_lines


    def _apply_hacks(self):
        if self.pdb_id:
            pdb_id = self.pdb_id.upper()
            if pdb_id == '2MBP':
                newlines = []
                added_COMPND = False
                for l in self.lines:
                    if l.startswith('COMPND'):
                        if not added_COMPND:
                            newlines.append('COMPND    MOL_ID: 1;')
                            newlines.append('COMPND   2 MOLECULE: MALTODEXTRIN-BINDING PROTEIN;')
                            newlines.append('COMPND   3 CHAIN: A;')
                            newlines.append('COMPND   4 ENGINEERED: YES')
                            added_COMPND = True
                    elif l.startswith("ATOM  ") or l.startswith("HETATM") or l.startswith("TER"):
                        newlines.append('%s%s%s' % (l[0:21], 'A', l[22:]))
                    elif l.startswith("SEQRES"):
                        newlines.append('%s%s%s' % (l[0:12], 'A', l[13:]))
                    else:
                        newlines.append(l)
                self.lines = newlines
            elif ROSETTA_HACKS_residues_to_remove.get(pdb_id):
                hacks = ROSETTA_HACKS_residues_to_remove[pdb_id]
                self.lines = [l for l in self.lines if not(l.startswith('ATOM'  )) or (l[21:27] not in hacks)]

        self._split_lines()

    ### Class functions ###

    @staticmethod
    def from_filepath(filepath, strict = True):
        '''A function to replace the old constructor call where a filename was passed in.'''
        return PDB(read_file(filepath), strict = strict)

    @staticmethod
    def from_lines(pdb_file_lines, strict = True):
        '''A function to replace the old constructor call where a list of the file's lines was passed in.'''
        return PDB("\n".join(pdb_file_lines), strict = strict)


    ### Private functions ###

    def _split_lines(self):
        '''Creates the parsed_lines dict which keeps all record data in document order indexed by the record type.'''
        parsed_lines = {}
        for rt in all_record_types:
            parsed_lines[rt] = []
        parsed_lines[0] = []

        for line in self.lines:
            linetype = line[0:6]
            if linetype in all_record_types:
                parsed_lines[linetype].append(line)
            else:
                parsed_lines[0].append(line)

        self.parsed_lines = parsed_lines
        self._update_structure_lines() # This does a second loop through the lines. We could do this logic above but I prefer a little performance hit for the cleaner code

    def _update_structure_lines(self):
        '''ATOM and HETATM lines may be altered by function calls. When this happens, this function should be called to keep self.structure_lines up to date.'''
        structure_lines = []
        atom_chain_order = []
        chain_atoms = {}

        for line in self.lines:
            linetype = line[0:6]
            if linetype == 'ATOM  ' or linetype == 'HETATM' or linetype == 'TER   ':
                chain_id = line[21]
                self.residue_types.add(line[17:20].strip())
                if missing_chain_ids.get(self.pdb_id):
                    chain_id = missing_chain_ids[self.pdb_id]
                structure_lines.append(line)
                if (chain_id not in atom_chain_order) and (chain_id != ' '):
                    atom_chain_order.append(chain_id)
                if linetype == 'ATOM  ':
                    atom_type = line[12:16].strip()
                    if atom_type:
                        chain_atoms[chain_id] = chain_atoms.get(chain_id, set())
                        chain_atoms[chain_id].add(atom_type)
            if linetype == 'ENDMDL':
                break

        self.structure_lines = structure_lines
        self.atom_chain_order = atom_chain_order
        self.chain_atoms = chain_atoms

    ### Basic functions ###

    def clone(self):
        '''A function to replace the old constructor call where a PDB object was passed in and 'cloned'.'''
        return PDB("\n".join(self.lines), pdb_id = self.pdb_id, strict = self.strict)

    def get_content(self):
        '''A function to replace the old constructor call where a PDB object was passed in and 'cloned'.'''
        return '\n'.join(self.lines)

    def write(self, pdbpath, separator = '\n'):
        write_file(pdbpath, separator.join(self.lines))

    def get_pdb_id(self):
        '''Return the PDB ID. If one was passed in to the constructor, this takes precedence, otherwise the header is
           parsed to try to find an ID. The header does not always contain a PDB ID in regular PDB files and appears to
           always have an ID of 'XXXX' in biological units so the constructor override is useful.'''
        if self.pdb_id:
            return self.pdb_id
        else:
            header = self.parsed_lines["HEADER"]
            assert(len(header) <= 1)
            if header:
                self.pdb_id = header[0][62:66]
                return self.pdb_id
        return None

    def get_ATOM_and_HETATM_chains(self):
        '''todo: remove this function as it now just returns a member element'''
        return self.atom_chain_order

    def get_annotated_chain_sequence_string(self, chain_id, use_seqres_sequences_if_possible, raise_Exception_if_not_found = True):
        '''A helper function to return the Sequence for a chain. If use_seqres_sequences_if_possible then we return the SEQRES
           Sequence if it exists. We return a tuple of values, the first identifying which sequence was returned.'''
        if use_seqres_sequences_if_possible and self.seqres_sequences and self.seqres_sequences.get(chain_id):
            return ('SEQRES', self.seqres_sequences[chain_id])
        elif self.atom_sequences.get(chain_id):
            return ('ATOM', self.atom_sequences[chain_id])
        elif raise_Exception_if_not_found:
            raise Exception('Error: Chain %s expected but not found.' % (str(chain_id)))
        else:
            return None

    def get_chain_sequence_string(self, chain_id, use_seqres_sequences_if_possible, raise_Exception_if_not_found = True):
        '''Similar to get_annotated_chain_sequence_string except that we only return the Sequence and do not state which sequence it was.'''
        chain_pair = self.get_annotated_chain_sequence_string(chain_id, use_seqres_sequences_if_possible, raise_Exception_if_not_found = raise_Exception_if_not_found)
        if chain_pair:
            return chain_pair[1]
        return None

    def _get_modified_residues(self):
        if not self.modified_residues:
            modified_residues = {}
            modified_residue_mapping_3 = {}

            # Add in the patch
            for k, v in modified_residues_patch.get(self.pdb_id, {}).iteritems():
                modified_residue_mapping_3[k] = v

            for line in self.parsed_lines["MODRES"]:
                restype = line[24:27].strip()
                restype_1 = residue_type_3to1_map.get(restype) or dna_nucleotides_2to1_map.get(restype)
                if not restype_1:
                    assert(restype in rna_nucleotides)
                    restype_1 = restype

                modified_residues["%s%s" % (line[16], line[18:23])] = {'modified_residue' : line[12:15], 'original_residue_3' : restype, 'original_residue_1' : restype_1}
                modified_residue_mapping_3[line[12:15]] = restype

            self.modified_residue_mapping_3 = modified_residue_mapping_3
            self.modified_residues = modified_residues


    def _get_replacement_pdb_id(self):
        '''Checks to see if the PDB file has been deprecated and, if so, what the new ID is.'''
        deprecation_lines = self.parsed_lines['OBSLTE']
        date_regex = re.compile('(\d+)-(\w{3})-(\d+)')
        if deprecation_lines:
            assert(len(deprecation_lines) == 1)
            tokens = deprecation_lines[0].split()[1:]
            if tokens[1].upper() in obsolete_pdb_ids_with_no_replacement_entries:
                assert(len(tokens) == 2)
            else:
                assert(len(tokens) == 3)
            if self.pdb_id:
                mtchs = date_regex.match(tokens[0])
                assert(mtchs)
                _day = int(mtchs.group(1))
                _month = mtchs.group(2)
                _year = int(mtchs.group(3)) # only two characters?
                assert(tokens[1] == self.pdb_id)
                self.deprecation_date = (_day, _month, _year) # no need to do anything fancier unless this is ever needed
                self.deprecated = True
                if len(tokens) == 3:
                    assert(len(tokens[2]) == 4)
                    self.replacement_pdb_id = tokens[2]

    ### PDB mutating functions ###

    def strip_to_chains(self, chains):
        '''Throw away all ATOM/HETATM/ANISOU/TER lines for chains that are not in the chains list.'''
        if chains:
            chains = set(chains)

            # Remove any structure lines associated with the chains
            self.lines = [l for l in self.lines if not(l.startswith('ATOM  ') or l.startswith('HETATM') or l.startswith('ANISOU') or l.startswith('TER')) or l[21] in chains]
            self._update_structure_lines()
            # todo: this logic should be fine if no other member elements rely on these lines e.g. residue mappings otherwise we need to update those elements here
        else:
            raise Exception('The chains argument needs to be supplied.')

    def strip_HETATMs(self, only_strip_these_chains = []):
        '''Throw away all HETATM lines. If only_strip_these_chains is specified then only strip HETATMs lines for those chains.'''
        if only_strip_these_chains:
            self.lines = [l for l in self.lines if not(l.startswith('HETATM')) or l[21] not in only_strip_these_chains]
        else:
            self.lines = [l for l in self.lines if not(l.startswith('HETATM'))]
        self._update_structure_lines()
        # todo: this logic should be fine if no other member elements rely on these lines e.g. residue mappings otherwise we need to update those elements here

    def generate_all_point_mutations_for_chain(self, chain_id):
        mutations = []
        if self.atom_sequences.get(chain_id):
            aas = sorted(residue_type_3to1_map.values())
            aas.remove('X')
            seq = self.atom_sequences[chain_id]
            for res_id in seq.order:
                r = seq.sequence[res_id]
                assert(chain_id == r.Chain)
                for mut_aa in aas:
                    if mut_aa != r.ResidueAA:
                        mutations.append(ChainMutation(r.ResidueAA, r.ResidueID, mut_aa, Chain = chain_id))
        return mutations

    ### FASTA functions ###


    def create_fasta(self, length = 80, prefer_seqres_order = True):
        fasta_string = ''
        if prefer_seqres_order:
            chain_order, sequences = self.seqres_chain_order or self.atom_chain_order, self.seqres_sequences or self.atom_sequences
        else:
            chain_order, sequences = self.atom_chain_order or self.seqres_chain_order, self.atom_sequences or self.seqres_sequences

        for c in chain_order:
            if c not in sequences:
                continue

            fasta_string += '>%s|%s|PDBID|CHAIN|SEQUENCE\n' % (self.pdb_id, c)
            seq = str(sequences[c])
            for line in [seq[x:x+length] for x in xrange(0, len(seq), length)]:
                fasta_string += line + '\n'

        return fasta_string


    ### PDB file parsing functions ###

    def _get_pdb_format_version(self):
        '''Remark 4 indicates the version of the PDB File Format used to generate the file.'''
        if not self.format_version:
            version = None
            version_lines = None
            try:
                version_lines = [line for line in self.parsed_lines['REMARK'] if int(line[7:10]) == 4 and line[10:].strip()]
            except: pass
            if version_lines:
                assert(len(version_lines) == 1)
                version_line = version_lines[0]
                version_regex = re.compile('.*?FORMAT V.(.*),')
                mtch = version_regex.match(version_line)
                if mtch and mtch.groups(0):
                    try:
                        version = float(mtch.groups(0)[0])
                    except:
                        pass
            self.format_version = version

    def get_resolution(self):
        resolution = None
        resolution_lines_exist = False
        for line in self.parsed_lines["REMARK"]:
            if line[9] == "2" and line[11:22] == "RESOLUTION.":
                #if id == :
                #	line = "REMARK   2 RESOLUTION. 3.00 ANGSTROMS.

                                # This code SHOULD work but there are badly formatted PDBs in the RCSB database.
                # e.g. "1GTX"
                #if line[31:41] == "ANGSTROMS.":
                #	try:
                #		resolution = float(line[23:30])
                #	except:
                #		raise Exception("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard. Expected data for diffraction experiments." % line )
                #if line[23:38] == "NOT APPLICABLE.":
                #	resolution = "N/A"
                #else:
                #	raise Exception("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard." % line )
                #
                # Instead, we use the code below:
                if resolution:
                    raise Exception("Found multiple RESOLUTION lines.")
                resolution_lines_exist = True
                strippedline = line[22:].strip()
                Aindex = strippedline.find("ANGSTROMS.")
                if strippedline == "NOT APPLICABLE.":
                    resolution = "N/A"
                elif Aindex != -1 and strippedline.endswith("ANGSTROMS."):
                    if strippedline[:Aindex].strip() == "NULL":
                        resolution = "N/A" # Yes, yes, yes, I know. Look at 1WSY.pdb.
                    else:
                        try:
                            resolution = float(strippedline[:Aindex].strip())
                        except:
                            raise PDBParsingException("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard. Expected data for diffraction experiments." % line )
                else:
                    raise PDBParsingException("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard." % line )
        if resolution_lines_exist and not resolution:
            raise PDBParsingException("Could not determine resolution.")
        return resolution

    def get_title(self):
        if self.parsed_lines.get("TITLE "):
            return " ".join([line[10:80].strip() for line in self.parsed_lines["TITLE "] if line[10:80].strip()])
        return None

    def get_techniques(self):
        techniques = None
        technique_lines_exist = False
        for line in self.parsed_lines["EXPDTA"]:
            technique_lines_exist = True
            techniques = line[10:71].split(";")
            for k in range(len(techniques)):
                techniques[k] = techniques[k].strip()
            techniques = ";".join(techniques)
        if technique_lines_exist and not techniques:
            raise PDBParsingException("Could not determine techniques used.")
        return techniques

    def get_UniProt_ACs(self):
        return [v['dbAccession'] for k, v in self.get_DB_references().get(self.pdb_id, {}).get('UNIPROT', {}).iteritems()]

    def get_DB_references(self):
        ''' "The DBREF record provides cross-reference links between PDB sequences (what appears in SEQRES record) and
                a corresponding database sequence." - http://www.wwpdb.org/documentation/format33/sect3.html#DBREF
        '''

        _database_names = {
            'GB'    :  'GenBank',
            'PDB'   :  'Protein Data Bank',
            'UNP'   :  'UNIPROT',
            'NORINE':  'Norine',
            'TREMBL': 'UNIPROT',
        }

        DBref = {}
        for l in self.parsed_lines["DBREF "]: # [l for l in self.lines if l.startswith('DBREF')]
            pdb_id = l[7:11]
            chain_id = l[12]
            seqBegin = int(l[14:18])
            insertBegin = l[18]
            seqEnd = int(l[20:24])
            insertEnd = l[24]
            database = _database_names[l[26:32].strip()]
            dbAccession = l[33:41].strip()
            dbIdCode = l[42:54].strip()
            dbseqBegin = int(l[55:60])
            idbnsBeg = l[60]
            dbseqEnd = int(l[62:67])
            dbinsEnd = l[67]

            DBref[pdb_id] = DBref.get(pdb_id, {})
            DBref[pdb_id][database] = DBref[pdb_id].get(database, {})
            if DBref[pdb_id][database].get(chain_id):
                if not(DBref[pdb_id][database][chain_id]['dbAccession'] == dbAccession and DBref[pdb_id][database][chain_id]['dbIdCode'] == dbIdCode):
                    raise PDBParsingException('This code needs to be generalized. dbIdCode should really be a list to handle chimera cases.')
            else:
                DBref[pdb_id][database][chain_id] = {'dbAccession'   :   dbAccession, 'dbIdCode'      :   dbIdCode, 'PDBtoDB_mapping' : []}

            DBref[pdb_id][database][chain_id]['PDBtoDB_mapping'].append(
                {'PDBRange'      :   ("%d%s" % (seqBegin,  insertBegin), "%d%s" % (seqEnd,  insertEnd)),
                'dbRange'       :   ("%d%s" % (dbseqBegin, idbnsBeg), "%d%s" % (dbseqEnd, dbinsEnd)),
                }
            )
        return DBref

    def get_molecules_and_source(self):
        # Check the COMPND lines
        COMPND_lines = self.parsed_lines["COMPND"]
        for x in range(1, len(COMPND_lines)):
            assert(int(COMPND_lines[x][7:10]) == x+1)
        if not COMPND_lines:
            raise MissingRecordsException("No COMPND records were found. Handle this gracefully.")

        # Concatenate the COMPND lines into one string, removing double spaces
        COMPND_lines = " ".join([line[10:].strip() for line in COMPND_lines])
        COMPND_lines.replace("  ", " ")

        # Split the COMPND lines into separate molecule entries
        molecules = {}
        MOL_DATA = ["MOL_ID:%s".strip() % s for s in COMPND_lines.split('MOL_ID:') if s]

        # Parse the molecule entries
        # The hacks below are due to some PDBs breaking the grammar by not following the standard which states:
        #   Specification: A String composed of a token and its associated value separated by a colon.
        #   Specification List: A sequence of Specifications, separated by semi-colons.
        # COMPND records are a specification list so semi-colons should not appear inside entries.
        # The hacks below could probably be removed if I assumed that the standard was not followed (valid) by
        #   e.g. splitting the COMPND data by allowed tokens (the keys of COMPND_field_map)
        # but I would want lots of tests in place first.
        for MD in MOL_DATA:
            # Hack for 2OMT
            MD = MD.replace('EPITHELIAL-CADHERIN; E-CAD/CTF1', 'EPITHELIAL-CADHERIN: E-CAD/CTF1')
            # Hack for 1M2T
            MD = MD.replace('SYNONYM: BETA-GALACTOSIDE SPECIFIC LECTIN I A CHAIN; MLA; ML-I A;', 'SYNONYM: BETA-GALACTOSIDE SPECIFIC LECTIN I A CHAIN, MLA, ML-I A,')
            # Hack for 1IBR
            MD = MD.replace('SYNONYM: RAN; TC4; RAN GTPASE; ANDROGEN RECEPTOR- ASSOCIATED PROTEIN 24;', 'SYNONYM: RAN TC4, RAN GTPASE, ANDROGEN RECEPTOR-ASSOCIATED PROTEIN 24;')
            # Hack for 1IBR
            MD = MD.replace('SYNONYM: KARYOPHERIN BETA-1 SUBUNIT; P95; NUCLEAR FACTOR P97; IMPORTIN 90', 'SYNONYM: KARYOPHERIN BETA-1 SUBUNIT, P95, NUCLEAR FACTOR P97, IMPORTIN 90')
            # Hack for 1NKH
            MD = MD.replace('SYNONYM: B4GAL-T1; BETA4GAL-T1; BETA-1,4-GALTASE 1; BETA-1, 4-GALACTOSYLTRANSFERASE 1;  UDP-GALACTOSE:BETA-N- ACETYLGLUCOSAMINE BETA-1,4-GALACTOSYLTRANSFERASE 1; EC: 2.4.1.22, 2.4.1.90, 2.4.1.38; ENGINEERED: YES; OTHER_DETAILS: CHAINS A AND B FORM FIRST, C AND D SECOND LACTOSE SYNTHASE COMPLEX',
                            'SYNONYM: B4GAL-T1, BETA4GAL-T1, BETA-1,4-GALTASE 1, BETA-1, 4-GALACTOSYLTRANSFERASE 1,  UDP-GALACTOSE:BETA-N- ACETYLGLUCOSAMINE BETA-1,4-GALACTOSYLTRANSFERASE 1, EC: 2.4.1.22, 2.4.1.90, 2.4.1.38, ENGINEERED: YES, OTHER_DETAILS: CHAINS A AND B FORM FIRST, C AND D SECOND LACTOSE SYNTHASE COMPLEX')
            # Hacks for 2PMI
            MD = MD.replace('SYNONYM: SERINE/THREONINE-PROTEIN KINASE PHO85; NEGATIVE REGULATOR OF THE PHO SYSTEM;',
                            'SYNONYM: SERINE/THREONINE-PROTEIN KINASE PHO85, NEGATIVE REGULATOR OF THE PHO SYSTEM;')
            MD = MD.replace('SYNONYM: PHOSPHATE SYSTEM CYCLIN PHO80; AMINOGLYCOSIDE ANTIBIOTIC SENSITIVITY PROTEIN 3;',
                            'SYNONYM: PHOSPHATE SYSTEM CYCLIN PHO80, AMINOGLYCOSIDE ANTIBIOTIC SENSITIVITY PROTEIN 3;')
            # Hack for 1JRH
            MD = MD.replace('FAB FRAGMENT;PEPSIN DIGESTION OF INTACT ANTIBODY', 'FAB FRAGMENT,PEPSIN DIGESTION OF INTACT ANTIBODY')
            # Hack for 1KJ1
            MD = MD.replace('SYNONYM: MANNOSE-SPECIFIC AGGLUTININ; LECGNA ', 'SYNONYM: MANNOSE-SPECIFIC AGGLUTININ, LECGNA ')
            # Hack for 1OCC - The Dean and I
            MD = MD.replace('SYNONYM: FERROCYTOCHROME C\:OXYGEN OXIDOREDUCTASE', 'SYNONYM: FERROCYTOCHROME C, OXYGEN OXIDOREDUCTASE')
            # Hack for 2AKY
            MD = MD.replace('SYNONYM: ATP\:AMP PHOSPHOTRANSFERASE, MYOKINASE', 'SYNONYM: ATP, AMP PHOSPHOTRANSFERASE, MYOKINASE')
            # Hack for 3BCI
            MD = MD.replace('SYNONYM: THIOL:DISULFIDE OXIDOREDUCTASE DSBA', 'SYNONYM: THIOL, DISULFIDE OXIDOREDUCTASE DSBA')
            # Hack for 3BCI
            MD = MD.replace('SYNONYM: THIOL:DISULFIDE OXIDOREDUCTASE DSBA', 'SYNONYM: THIOL, DISULFIDE OXIDOREDUCTASE DSBA')
            # Hack for 1ELV
            MD = MD.replace('FRAGMENT: CCP2-SP CATALYTIC FRAGMENT: ASP363-ASP-673 SEGMENT PRECEDED BY AN ASP-LEU SEQUENCE ADDED AT THE N-TERMINAL END',
                            'FRAGMENT: CCP2-SP CATALYTIC FRAGMENT; ASP363-ASP-673 SEGMENT PRECEDED BY AN ASP-LEU SEQUENCE ADDED AT THE N-TERMINAL END')
            # Hack for 1E6E
            MD = MD.replace('MOLECULE: NADPH\:ADRENODOXIN OXIDOREDUCTASE;', 'MOLECULE: NADPH;ADRENODOXIN OXIDOREDUCTASE;')
            # Hack for 1JZD
            MD = MD.replace('MOLECULE: THIOL:DISULFIDE INTERCHANGE PROTEIN', 'MOLECULE: THIOL;DISULFIDE INTERCHANGE PROTEIN')
            # Hack for 1N2C
            MD = MD.replace('OTHER_DETAILS: 2\:1 COMPLEX OF HOMODIMERIC FE-PROTEIN', 'OTHER_DETAILS: 2;1 COMPLEX OF HOMODIMERIC FE-PROTEIN')
            # Hack for 1S6P
            MD = MD.replace('MOLECULE: POL POLYPROTEIN [CONTAINS: REVERSE TRANSCRIPTASE]', 'MOLECULE: POL POLYPROTEIN [CONTAINS; REVERSE TRANSCRIPTASE]')
            # Hack for 1Z9E
            MD = MD.replace('FRAGMENT: SEQUENCE DATABASE RESIDUES 347-471 CONTAINS: HIV- 1 INTEGRASE-BINDING DOMAIN', 'FRAGMENT: SEQUENCE DATABASE RESIDUES 347-471 CONTAINS; HIV- 1 INTEGRASE-BINDING DOMAIN')
            # Hacks for 2GOX
            MD = MD.replace('FRAGMENT: FRAGMENT OF ALPHA CHAIN: RESIDUES 996-1287;', 'FRAGMENT: FRAGMENT OF ALPHA CHAIN; RESIDUES 996-1287;')
            MD = MD.replace('FRAGMENT: C-TERMINAL DOMAIN: RESIDUES 101-165;', 'FRAGMENT: C-TERMINAL DOMAIN; RESIDUES 101-165;')

            MOL_fields = [s.strip() for s in MD.split(';') if s.strip()]

            molecule = {}
            for field in MOL_fields:
                field = field.split(":")
                assert(1 <= len(field) <= 2)
                if len(field) == 2: # Hack for 1MBG - missing field value
                    field_name = COMPND_field_map[field[0].strip()]
                    field_data = field[1].strip()
                    molecule[field_name] = field_data

            ### Normalize and type the fields ###

            # Required (by us) fields
            molecule['MoleculeID'] = int(molecule['MoleculeID'])
            molecule['Chains'] = map(string.strip, molecule['Chains'].split(','))
            for c in molecule['Chains']:
                assert(len(c) == 1)

            # Optional fields
            if not molecule.get('Engineered'):
                molecule['Engineered'] = None
            elif molecule.get('Engineered') == 'YES':
                molecule['Engineered'] = True
            elif molecule.get('Engineered') == 'NO':
                molecule['Engineered'] = False
            else:
                raise PDBParsingException("Error parsing ENGINEERED field of COMPND lines. Expected 'YES' or 'NO', got '%s'." % molecule['Engineered'])

            if molecule.get('Mutation'):
                if molecule['Mutation'] != 'YES':
                    raise PDBParsingException("Error parsing MUTATION field of COMPND lines. Expected 'YES', got '%s'." % molecule['Mutation'])
                else:
                    molecule['Mutation'] = True
            else:
                molecule['Mutation'] = None

            # Add missing fields
            for k in COMPND_field_map.values():
                if k not in molecule.keys():
                    molecule[k] = None

            molecules[molecule['MoleculeID']] = molecule

        # Extract the SOURCE lines
        SOURCE_lines = self.parsed_lines["SOURCE"]
        for x in range(1, len(SOURCE_lines)):
            assert(int(SOURCE_lines[x][7:10]) == x+1)
        if not SOURCE_lines:
            raise MissingRecordsException("No SOURCE records were found. Handle this gracefully.")

        # Concatenate the SOURCE lines into one string, removing double spaces
        SOURCE_lines = " ".join([line[10:].strip() for line in SOURCE_lines])
        SOURCE_lines.replace("  ", " ")

        # Split the SOURCE lines into separate molecule entries
        MOL_DATA = ["MOL_ID:%s".strip() % s for s in SOURCE_lines.split('MOL_ID:') if s]
        # Parse the molecule entries
        for MD in MOL_DATA:
            MOL_fields = [s.strip() for s in MD.split(';') if s.strip()]
            new_molecule = {}
            for field in MOL_fields:
                field = field.split(":")
                if SOURCE_field_map.get(field[0].strip()):
                    field_name = SOURCE_field_map[field[0].strip()]
                    field_data = field[1].strip()
                    new_molecule[field_name] = field_data

            MoleculeID = int(new_molecule['MoleculeID'])
            assert(MoleculeID in molecules)
            molecule = molecules[MoleculeID]

            for field_name, field_data in new_molecule.iteritems():
                if field_name != 'MoleculeID':
                    molecule[field_name] = field_data

            # Normalize and type the fields

            if not molecule.get('Synthetic'):
                molecule['Synthetic'] = None
            elif molecule.get('Synthetic') == 'YES':
                molecule['Synthetic'] = True
            elif molecule.get('Synthetic') == 'NO':
                molecule['Synthetic'] = False
            else:
                raise PDBParsingException("Error parsing SYNTHETIC field of SOURCE lines. Expected 'YES' or 'NO', got '%s'." % molecule['Synthetic'])

            # Add missing fields
            for k in SOURCE_field_map.values():
                if k not in molecule.keys():
                    molecule[k] = None

        return [v for k, v in sorted(molecules.iteritems())]

    def get_journal(self):
        if self.parsed_lines["JRNL  "]:
            if not self.journal:
                self.journal = JRNL(self.parsed_lines["JRNL  "])
            return self.journal.get_info()
        return None

    ### Sequence-related functions ###

    def _get_SEQRES_sequences(self):
        '''Creates the SEQRES Sequences and stores the chains in order of their appearance in the SEQRES records. This order of chains
           in the SEQRES sequences does not always agree with the order in the ATOM records.'''

        pdb_id = self.get_pdb_id()
        SEQRES_lines = self.parsed_lines["SEQRES"]

        modified_residue_mapping_3 = self.modified_residue_mapping_3
        # I commented this out since we do not need it for my current test cases
        #for k, v in self.modified_residues.iteritems():
        #    assert(v['modified_residue'] not in modified_residues)
        #    modified_residues[v['modified_residue']] = v['original_residue_3']

        for x in range(0, len(SEQRES_lines)):
            assert(SEQRES_lines[x][7:10].strip().isdigit())

        seqres_chain_order = []
        SEQRES_lines = [line[11:].rstrip() for line in SEQRES_lines] # we cannot strip the left characters as some cases e.g. 2MBP are missing chain identifiers

        # Collect all residues for all chains, remembering the chain order
        chain_tokens = {}
        for line in SEQRES_lines:
            chainID = line[0]
            if missing_chain_ids.get(self.pdb_id):
                chainID = missing_chain_ids[self.pdb_id]
            if chainID not in seqres_chain_order:
                seqres_chain_order.append(chainID)
            chain_tokens[chainID] = chain_tokens.get(chainID, [])
            chain_tokens[chainID].extend(line[6:].strip().split())

        sequences = {}
        self.chain_types = {}

        for chain_id, tokens in chain_tokens.iteritems():

            # Determine whether the chain is DNA, RNA, or a protein chain
            # 1H38 is a good test for this - it contains DNA (chains E and G and repeated by H, K, N, J, M, P), RNA (chain F, repeated by I, L, O) and protein (chain D, repeated by A,B,C) sequences
            # 1ZC8 is similar but also has examples of DU
            # 4IHY has examples of DI (I is inosine)
            # 2GRB has RNA examples of I and U
            # 1LRP has protein chains with only CA atoms
            # This will throw an exception when a non-canonical is found which is not listed in basics.py. In that case, the list in basics.py should be updated.

            chain_type = None
            set_of_tokens = set(tokens)
            if (set(tokens).union(all_recognized_dna) == all_recognized_dna):# or (len(set_of_tokens) <= 5 and len(set_of_tokens.union(dna_nucleotides)) == len(set_of_tokens) + 1): # allow one unknown DNA residue
                chain_type = 'DNA'
            elif (set(tokens).union(all_recognized_rna) == all_recognized_rna):# or (len(set_of_tokens) <= 5 and len(set_of_tokens.union(dna_nucleotides)) == len(set_of_tokens) + 1): # allow one unknown DNA residue
                chain_type = 'RNA'
            else:
                assert(len(set(tokens).intersection(dna_nucleotides)) == 0)
                assert(len(set(tokens).intersection(rna_nucleotides)) == 0)
                chain_type = 'Protein'
                if not self.chain_atoms.get(chain_id):
                    # possible for biological unit files
                    continue
                if self.chain_atoms[chain_id] == set(['CA']):
                    chain_type = 'Protein skeleton'

            # Get the sequence, mapping non-canonicals to the appropriate letter
            self.chain_types[chain_id] = chain_type
            sequence = []
            if chain_type == 'DNA':
                for r in tokens:
                    if dna_nucleotides_2to1_map.get(r):
                        sequence.append(dna_nucleotides_2to1_map[r])
                    else:
                        if non_canonical_dna.get(r):
                            sequence.append(non_canonical_dna[r])
                        else:
                            raise Exception("Unknown DNA residue %s." % r)
            elif chain_type == 'RNA':
                for r in tokens:
                    if r in rna_nucleotides:
                        sequence.append(r)
                    else:
                        if non_canonical_rna.get(r):
                            sequence.append(non_canonical_rna[r])
                        else:
                            raise Exception("Unknown RNA residue %s." % r)
            else:
                token_counter = 0
                for r in tokens:
                    token_counter += 1
                    if residue_type_3to1_map.get(r):
                        sequence.append(residue_type_3to1_map[r])
                    else:

                        if self.modified_residue_mapping_3.get(r):
                            sequence.append(residue_type_3to1_map[self.modified_residue_mapping_3.get(r)])
                        elif non_canonical_amino_acids.get(r):
                            #print('Mapping non-canonical residue %s to %s.' % (r, non_canonical_amino_acids[r]))
                            #print(SEQRES_lines)
                            #print(line)
                            sequence.append(non_canonical_amino_acids[r])
                        elif r == 'UNK':
                            continue
                        # Skip these residues
                        elif r == 'ACE' and token_counter == 1:
                            # Always allow ACE as the first residue of a chain
                            sequence.append('X')
                        elif r == 'ACE' and pdb_id in cases_with_ACE_residues_we_can_ignore:
                            sequence.append('X')
                            #continue
                        # End of skipped residues
                        else:
                            #print(modified_residue_mapping_3)
                            if modified_residue_mapping_3.get(r):
                                if modified_residue_mapping_3[r] == 'UNK':
                                    sequence.append('X')
                                else:
                                    assert(modified_residue_mapping_3[r] in residue_types_3)
                                    sequence.append(residue_type_3to1_map[modified_residue_mapping_3[r]])
                            else:
                                raise Exception("Unknown protein residue %s in chain %s." % (r, chain_id))
            sequences[chain_id] = "".join(sequence)

        self.seqres_chain_order = seqres_chain_order

        # Create Sequence objects for the SEQRES sequences
        for chain_id, sequence in sequences.iteritems():
            self.seqres_sequences[chain_id] = Sequence.from_sequence(chain_id, sequence, self.chain_types[chain_id])


    def _get_ATOM_sequences(self):
        '''Creates the ATOM Sequences.'''

        # Get a list of all residues with ATOM or HETATM records
        atom_sequences = {}
        structural_residue_IDs_set = set() # use a set for a quicker lookup
        ignore_HETATMs = True # todo: fix this if we need to deal with HETATMs

        residue_lines_by_chain = []
        structural_residue_IDs_set = []

        model_index = 0
        residue_lines_by_chain.append([])
        structural_residue_IDs_set.append(set())
        full_code_map = {}
        for l in self.structure_lines:
            if l.startswith("TER   "):
                model_index += 1
                residue_lines_by_chain.append([])
                structural_residue_IDs_set.append(set())
            else:
                residue_id = l[21:27]
                if residue_id not in structural_residue_IDs_set[model_index]:
                    residue_lines_by_chain[model_index].append(l)
                    structural_residue_IDs_set[model_index].add(residue_id)
                full_code_map[l[21]] = full_code_map.get(l[21], set())
                full_code_map[l[21]].add(l[17:20].strip())

        # Get the residues used by the residue lines. These can be used to determine the chain type if the header is missing.
        for chain_id in self.atom_chain_order:
            if full_code_map.get(chain_id):
                # The chains may contain other molecules e.g. MG or HOH so before we decide their type based on residue types alone,
                # we subtract out those non-canonicals
                canonical_molecules = full_code_map[chain_id].intersection(dna_nucleotides.union(rna_nucleotides).union(residue_types_3))
                if canonical_molecules.union(dna_nucleotides) == dna_nucleotides:
                    self.chain_types[chain_id] = 'DNA'
                elif canonical_molecules.union(rna_nucleotides) == rna_nucleotides:
                    self.chain_types[chain_id] = 'RNA'
                else:
                    self.chain_types[chain_id] = 'Protein'

        line_types_by_chain = []
        chain_ids = []
        for model_index in range(len(residue_lines_by_chain)):
            line_types = set()
            if residue_lines_by_chain[model_index]:
                if missing_chain_ids.get(self.pdb_id):
                    chain_ids.append(missing_chain_ids[self.pdb_id])
                else:
                    chain_ids.append(residue_lines_by_chain[model_index][0][21])
            for l in residue_lines_by_chain[model_index]:
                line_types.add(l[0:6])
            if line_types == set(['ATOM']):
                line_types_by_chain.append('ATOM')
            elif line_types == set(['HETATM']):
                line_types_by_chain.append('HETATM')
            else:
                line_types_by_chain.append('Mixed')

        for x in range(0, len(residue_lines_by_chain)):
            residue_lines = residue_lines_by_chain[x]
            line_types = line_types_by_chain[x]
            if ignore_HETATMs and line_types == 'HETATM':
                continue

            for y in range(len(residue_lines)):
                l = residue_lines[y]
                residue_type = l[17:20].strip()
                if l.startswith("HETATM"):
                    if self.modified_residue_mapping_3.get(residue_type):
                        residue_type = self.modified_residue_mapping_3[residue_type]
                    elif y == (len(residue_lines) - 1):
                        # last residue in the chain
                        if residue_type == 'NH2':
                            residue_type = 'UNK' # fixes a few cases e.g. 1MBG, 1K9Q, 1KA6
                        elif ignore_HETATMs:
                            continue

                    elif ignore_HETATMs:
                        continue

                residue_id = l[21:27]
                chain_id = l[21]
                if missing_chain_ids.get(self.pdb_id):
                    chain_id = missing_chain_ids[self.pdb_id]

                if chain_id in self.chain_types:
                    # This means the pdb had SEQRES and we constructed atom_sequences
                    chain_type = self.chain_types[chain_id]
                else:
                    # Otherwise assume this is protein
                    chain_type = 'Protein'

                atom_sequences[chain_id] = atom_sequences.get(chain_id, Sequence(chain_type))

                residue_type = self.modified_residue_mapping_3.get(residue_type, residue_type)

                short_residue_type = None
                if residue_type == 'UNK':
                    short_residue_type = 'X'
                elif chain_type == 'Protein' or chain_type == 'Protein skeleton':
                    short_residue_type = residue_type_3to1_map.get(residue_type) or protonated_residue_type_3to1_map.get(residue_type) or non_canonical_amino_acids.get(residue_type)
                elif chain_type == 'DNA':
                    short_residue_type = dna_nucleotides_2to1_map.get(residue_type) or non_canonical_dna.get(residue_type)
                elif chain_type == 'RNA':
                    short_residue_type = non_canonical_rna.get(residue_type) or residue_type

                if not short_residue_type:
                    if l.startswith("ATOM") and l[12:16] == ' OH2' and l[17:20] == 'TIP':
                        continue
                    elif not self.strict:
                        short_residue_type = 'X'
                    else:
                        raise NonCanonicalResidueException("Unrecognized residue type %s in PDB file '%s', residue ID '%s'." % (residue_type, str(self.pdb_id), str(residue_id)))

                #structural_residue_IDs.append((residue_id, short_residue_type))
                # KAB - way to allow for multiresidue noncanonical AA's
                if len(short_residue_type) == 1:
                    atom_sequences[chain_id].add(PDBResidue(residue_id[0], residue_id[1:], short_residue_type, chain_type))
                else:
                    for char in short_residue_type:
                        atom_sequences[chain_id].add(PDBResidue(residue_id[0], residue_id[1:], char, chain_type))

        self.atom_sequences = atom_sequences

    def _get_ATOM_sequences_2(self):
        '''Creates the ATOM Sequences.'''

        # Get a list of all residues with ATOM or HETATM records
        atom_sequences = {}
        structural_residue_IDs_set = set() # use a set for a quicker lookup
        ignore_HETATMs = True # todo: fix this if we need to deal with HETATMs
        for l in self.structure_lines:
            residue_type = l[17:20].strip()
            if l.startswith("HETATM"):
                if self.modified_residue_mapping_3.get(residue_type):
                    residue_type = self.modified_residue_mapping_3[residue_type]
                elif ignore_HETATMs:
                    continue

            residue_id = l[21:27]
            if residue_id not in structural_residue_IDs_set:
                chain_id = l[21]
                chain_type = self.chain_types[chain_id]
                atom_sequences[chain_id] = atom_sequences.get(chain_id, Sequence(chain_type))
                residue_type = l[17:20].strip()

                residue_type = self.modified_residue_mapping_3.get(residue_type, residue_type)
                short_residue_type = None
                if residue_type == 'UNK':
                    short_residue_type = 'X'
                elif chain_type == 'Protein' or chain_type == 'Protein skeleton':
                    short_residue_type = residue_type_3to1_map.get(residue_type) or protonated_residue_type_3to1_map.get(residue_type)
                elif chain_type == 'DNA':
                    short_residue_type = dna_nucleotides_2to1_map.get(residue_type) or non_canonical_dna.get(residue_type)
                elif chain_type == 'RNA':
                    short_residue_type = non_canonical_rna.get(residue_type) or residue_type
                elif not self.strict:
                    short_residue_type = 'X'
                else:
                    raise NonCanonicalResidueException("Unrecognized residue type %s in PDB file '%s', residue ID '%s'." % (residue_type, str(self.pdb_id), str(residue_id)))

                #structural_residue_IDs.append((residue_id, short_residue_type))
                atom_sequences[chain_id].add(PDBResidue(residue_id[0], residue_id[1:], short_residue_type, chain_type))
                structural_residue_IDs_set.add(residue_id)

        self.atom_sequences = atom_sequences


    def construct_pdb_to_rosetta_residue_map(self, rosetta_scripts_path, rosetta_database_path, extra_command_flags = None):
        ''' Uses the features database to create a mapping from Rosetta-numbered residues to PDB ATOM residues.
            Next, the object's rosetta_sequences (a dict of Sequences) element is created.
            Finally, a SequenceMap object is created mapping the Rosetta Sequences to the ATOM Sequences.

            The extra_command_flags parameter expects a string e.g. "-ignore_zero_occupancy false".
        '''

        ## Create a mapping from Rosetta-numbered residues to PDB ATOM residues

        # Apply any PDB-specific hacks
        specific_flag_hacks = None
        if self.pdb_id and HACKS_pdb_specific_hacks.get(self.pdb_id):
            specific_flag_hacks = HACKS_pdb_specific_hacks[self.pdb_id]

        skeletal_chains = sorted([k for k in self.chain_types.keys() if self.chain_types[k] == 'Protein skeleton'])
        if skeletal_chains:
            raise PDBMissingMainchainAtomsException('The PDB to Rosetta residue map could not be created as chains %s only have CA atoms present.' % ", ".join(skeletal_chains))

        # Get the residue mapping using the features database
        pdb_file_contents = "\n".join(self.structure_lines)
        success, mapping = get_pdb_contents_to_pose_residue_map(pdb_file_contents, rosetta_scripts_path, rosetta_database_path, pdb_id = self.pdb_id, extra_flags = ((specific_flag_hacks or '') + ' ' + (extra_command_flags or '')).strip())
        if not success:
            raise Exception("An error occurred mapping the PDB ATOM residue IDs to the Rosetta numbering.\n%s" % "\n".join(mapping))

        ## Create Sequences for the Rosetta residues (self.rosetta_sequences)

        # Initialize maps
        rosetta_residues = {}
        rosetta_sequences = {}
        for chain_id in self.atom_chain_order:
            chain_type = self.chain_types[chain_id]
            rosetta_residues[chain_id] = {}
            rosetta_sequences[chain_id] = Sequence(chain_type)

        # Create a map rosetta_residues, Chain -> Rosetta residue ID -> Rosetta residue information
        rosetta_pdb_mappings = {}
        for chain_id in self.atom_chain_order:
            rosetta_pdb_mappings[chain_id] = {}
        for k, v in mapping.iteritems():
            rosetta_residues[k[0]][v['pose_residue_id']] = v
            rosetta_pdb_mappings[k[0]][v['pose_residue_id']] = k

        # Create rosetta_sequences map Chain -> Sequence(Residue)
        for chain_id, v in sorted(rosetta_residues.iteritems()):
            chain_type = self.chain_types[chain_id]
            for rosetta_id, residue_info in sorted(v.iteritems()):
                short_residue_type = None

                if chain_type == 'Protein':
                    residue_type = residue_info['name3'].strip()
                    short_residue_type = residue_type_3to1_map[residue_type]
                else:
                    assert(chain_type == 'DNA' or chain_type == 'RNA')
                    residue_type = residue_info['res_type'].strip()
                    if residue_type.find('UpperDNA') != -1 or residue_type.find('LowerDNA') != -1:
                        residue_type = residue_type[:3]
                    short_residue_type = dna_nucleotides_3to1_map.get(residue_type) # Commenting this out since Rosetta does not seem to handle these "or non_canonical_dna.get(residue_type)"

                assert(short_residue_type)
                rosetta_sequences[chain_id].add(Residue(chain_id, rosetta_id, short_residue_type, chain_type))


        ## Create SequenceMap objects to map the Rosetta Sequences to the ATOM Sequences
        rosetta_to_atom_sequence_maps = {}
        for chain_id, rosetta_pdb_mapping in rosetta_pdb_mappings.iteritems():
            rosetta_to_atom_sequence_maps[chain_id] = SequenceMap.from_dict(rosetta_pdb_mapping)

        self.rosetta_to_atom_sequence_maps = rosetta_to_atom_sequence_maps
        self.rosetta_sequences = rosetta_sequences


    def get_atom_sequence_to_rosetta_map(self):
        '''Uses the Rosetta->ATOM injective map to construct an injective mapping from ATOM->Rosetta.
           We do not extend the injection to include ATOM residues which have no corresponding Rosetta residue.
             e.g. atom_sequence_to_rosetta_mapping[c].map.get('A  45 ') will return None if there is no corresponding Rosetta residue
           those residues to None.
           Likewise, if a PDB chain c is not present in the Rosetta model then atom_sequence_to_rosetta_mapping[c].map.get(s) returns None.
        '''
        if not self.rosetta_to_atom_sequence_maps and self.rosetta_sequences:
            raise Exception('The PDB to Rosetta mapping has not been determined. Please call construct_pdb_to_rosetta_residue_map first.')

        atom_sequence_to_rosetta_mapping = {}
        for chain_id, mapping in self.rosetta_to_atom_sequence_maps.iteritems():
            chain_mapping = {}
            for k in mapping:
                chain_mapping[k[1]] = k[0]
            atom_sequence_to_rosetta_mapping[chain_id] = SequenceMap.from_dict(chain_mapping)

        # Add empty maps for missing chains
        for chain_id, sequence in self.atom_sequences.iteritems():
            if not atom_sequence_to_rosetta_mapping.get(chain_id):
                atom_sequence_to_rosetta_mapping[chain_id] = SequenceMap()

        return atom_sequence_to_rosetta_mapping


    def get_atom_sequence_to_rosetta_json_map(self):
        '''Returns the mapping from PDB ATOM residue IDs to Rosetta residue IDs in JSON format.'''
        import json
        d = {}
        atom_sequence_to_rosetta_mapping = self.get_atom_sequence_to_rosetta_map()
        for c, sm in atom_sequence_to_rosetta_mapping.iteritems():
            for k, v in sm.map.iteritems():
                d[k] = v
        return json.dumps(d, sort_keys = True)


    def get_rosetta_sequence_to_atom_json_map(self):
        '''Returns the mapping from Rosetta residue IDs to PDB ATOM residue IDs in JSON format.'''
        import json
        if not self.rosetta_to_atom_sequence_maps and self.rosetta_sequences:
            raise Exception('The PDB to Rosetta mapping has not been determined. Please call construct_pdb_to_rosetta_residue_map first.')

        d = {}
        for c, sm in self.rosetta_to_atom_sequence_maps.iteritems():
            for k, v in sm.map.iteritems():
                d[k] = v
            #d[c] = sm.map
        return json.dumps(d, sort_keys = True)


    def map_pdb_residues_to_rosetta_residues(self, mutations):
        '''This function takes a list of ChainMutation objects and uses the PDB to Rosetta mapping to return the corresponding
           list of SimpleMutation objects using Rosetta numbering.
           e.g.
              p = PDB(...)
              p.construct_pdb_to_rosetta_residue_map()
              rosetta_mutations = p.map_pdb_residues_to_rosetta_residues(pdb_mutations)
        '''
        if not self.rosetta_to_atom_sequence_maps and self.rosetta_sequences:
            raise Exception('The PDB to Rosetta mapping has not been determined. Please call construct_pdb_to_rosetta_residue_map first.')

        rosetta_mutations = []
        atom_sequence_to_rosetta_mapping = self.get_atom_sequence_to_rosetta_map()
        for m in mutations:
            rosetta_residue_id = atom_sequence_to_rosetta_mapping[m.Chain].get('%s%s' % (m.Chain, m.ResidueID))
            rosetta_mutations.append(SimpleMutation(m.WildTypeAA, rosetta_residue_id, m.MutantAA))
        return rosetta_mutations


    def assert_wildtype_matches(self, mutation):
        '''Check that the wildtype of the Mutation object matches the PDB sequence.'''
        readwt = self.getAminoAcid(self.getAtomLine(mutation.Chain, mutation.ResidueID))
        assert(mutation.WildTypeAA == residue_type_3to1_map[readwt])


    ### END OF REFACTORED CODE


    @staticmethod
    def ChainResidueID2String(chain, residueID):
        '''Takes a chain ID e.g. 'A' and a residueID e.g. '123' or '123A' and returns the 6-character identifier spaced as in the PDB format.'''
        return "%s%s" % (chain, PDB.ResidueID2String(residueID))

    @staticmethod
    def ResidueID2String(residueID):
        '''Takes a residueID e.g. '123' or '123A' and returns the 5-character identifier spaced as in the PDB format.'''
        if residueID.isdigit():
            return "%s " % (residueID.rjust(4))
        else:
            return "%s" % (residueID.rjust(5))

    def validate_mutations(self, mutations):
        '''This function has been refactored to use the SimpleMutation class.
           The parameter is a list of Mutation objects. The function has no return value but raises a PDBValidationException
           if the wildtype in the Mutation m does not match the residue type corresponding to residue m.ResidueID in the PDB file.
        '''
        # Chain, ResidueID, WildTypeAA, MutantAA
        resID2AA = self.get_residue_id_to_type_map()
        badmutations = []
        for m in mutations:
            wildtype = resID2AA.get(PDB.ChainResidueID2String(m.Chain, m.ResidueID), "")
            if m.WildTypeAA != wildtype:
                badmutations.append(m)
        if badmutations:
            raise PDBValidationException("The mutation(s) %s could not be matched against the PDB %s." % (", ".join(map(str, badmutations)), self.pdb_id))


    def getAminoAcid(self, line):
        return line[17:20]


    def getAtomLine(self, chain, resid):
        '''This function assumes that all lines are ATOM or HETATM lines.
           resid should have the proper PDB format i.e. an integer left-padded
           to length 4 followed by the insertion code which may be a blank space.'''
        for line in self.lines:
            fieldtype = line[0:6].strip()
            assert(fieldtype == "ATOM" or fieldtype == "HETATM")
            if line[21:22] == chain and resid == line[22:27]:
                return line
        raise Exception("Could not find the ATOM/HETATM line corresponding to chain '%(chain)s' and residue '%(resid)s'." % vars())


    def CheckForPresenceOf(self, reslist):
        '''This checks whether residues in reslist exist in the ATOM lines.
           It returns a list of the residues in reslist which did exist.'''
        if type(reslist) == type(""):
            reslist = [reslist]

        foundRes = {}
        for line in self.lines:
            resname = line[17:20]
            if line[0:4] == "ATOM":
                if resname in reslist:
                    foundRes[resname] = True

        return foundRes.keys()


    def get_residue_id_to_type_map(self):
        '''Returns a dictionary mapping 6-character residue IDs (Chain, residue number, insertion code e.g. "A 123B") to the
           corresponding one-letter amino acid.

           Caveat: This function ignores occupancy - this function should be called once occupancy has been dealt with appropriately.'''

        resid2type = {}
        atomlines = self.parsed_lines['ATOM  ']
        for line in atomlines:
            resname = line[17:20]
            if resname in allowed_PDB_residues_types and line[13:16] == 'CA ':
                resid2type[line[21:27]] = residue_type_3to1_map.get(resname) or protonated_residue_type_3to1_map.get(resname)
        return resid2type


    def chain_ids(self):
        chain_ids = set()
        chainlist = []
        for line in self.lines:
            if line[0:4] == "ATOM" and line[17:20] in allowed_PDB_residues_types and line[26] == ' ':
                chain = line[21:22]
                if chain not in chain_ids:
                    chain_ids.add(chain)
                    chainlist.append(chain)

        return chainlist
