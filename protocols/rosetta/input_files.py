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
from basics import SimpleMutation


class RosettaFileParsingException(Exception): pass


class Mutfile (object):
    '''Note: This class behaves differently to Resfile. It stores mutation information using the SimpleMutation class.

       Rosetta mutfiles are text files split into sections where each section contains a number of mutations i.e. each
       section defines one mutagenesis.
       Mutfile objects represent the contents of these files by storing the mutations as a list of SimpleMutation lists.
    '''

    header_pattern = '^total\s+(\d+)\s*(?:#.*)?$'
    mutation_group_header_pattern = '^\s*(\d+)\s*(?:#.*)?$'
    mutation_pattern = '^\s*([A-Z])\s+(\d+)\s+([A-Z])\s*(?:#.*)?$'


    @staticmethod
    def from_file(filepath):
        return Mutfile(open(filepath).read())


    @staticmethod
    def from_mutagenesis(mutations):
        '''This is a special case (the common case) of from_mutations where there is only one mutagenesis/mutation group.'''
        return Mutfile.from_mutageneses([mutations])


    @staticmethod
    def from_mutageneses(mutation_groups):
        '''mutation_groups is expected to be a list containing lists of SimpleMutation objects.'''
        mf = Mutfile()
        mf.mutation_groups = mutation_groups
        return mf


    def __repr__(self):
        '''Creates a mutfile from the set of mutation groups.'''
        s = []

        # Header
        total_number_of_mutations = sum([len(mg) for mg in self.mutation_groups])
        s.append('total %d' % total_number_of_mutations)

        # Mutation groups
        for mg in self.mutation_groups:
            assert(len(mg) > 0)
            s.append('%d' % len(mg))
            # Mutation list
            for m in mg:
                s.append('%(WildTypeAA)s %(ResidueID)d %(MutantAA)s' % m.__dict__)

        s.append('')
        return '\n'.join(s)


    def __init__(self, mutfile_content = None):

        self.mutation_groups = []
        if mutfile_content:
            # Parse the file header
            mutfile_content = mutfile_content.strip()
            data_lines = [l for l in mutfile_content.split('\n') if l.strip()]
            try:
                num_mutations = int(re.match(Mutfile.header_pattern, data_lines[0]).group(1))
            except:
                raise RosettaFileParsingException('The mutfile has a bad header (expected "total n" where n is an integer).')

            line_counter, mutation_groups = 1, []
            while True:
                if line_counter >= len(data_lines):
                    break

                mutation_group_number = len(mutation_groups) + 1

                # Parse the group header
                try:
                    group_header = data_lines[line_counter]
                    line_counter += 1
                    num_mutations_in_group = int(re.match(Mutfile.mutation_group_header_pattern, group_header).group(1))
                    if num_mutations_in_group < 1:
                        raise RosettaFileParsingException('The mutfile has a record in mutation group %d: the number of reported mutations must be an integer greater than zero.' % mutation_group_number)
                except:
                    raise RosettaFileParsingException('The mutfile has a bad header for mutation group %d.' % mutation_group_number)

                # Parse the mutations in the group
                try:
                    mutations = []
                    for mutation_line in data_lines[line_counter: line_counter + num_mutations_in_group]:
                        mtch = re.match(Mutfile.mutation_pattern, mutation_line)
                        mutations.append(SimpleMutation(mtch.group(1), int(mtch.group(2)), mtch.group(3)))
                    mutation_groups.append(mutations)
                    line_counter += num_mutations_in_group
                except:
                    raise RosettaFileParsingException('An exception occurred while parsing the mutations for mutation group %d.' % mutation_group_number)

            if sum([len(mg) for mg in mutation_groups]) != num_mutations:
                raise RosettaFileParsingException('A total of %d mutations were expected from the file header but the file contained %d mutations.' % (num_mutations, sum([len(mg) for mg in mutation_groups])))

            self.mutation_groups = mutation_groups


    def get_total_mutation_count(self):
        return sum([len(mg) for mg in self.mutation_groups])
