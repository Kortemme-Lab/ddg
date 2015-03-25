# The MIT License (MIT)
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

kortemme_baker_protein_protein_interfaces = [
#1A22
#Kortemme & Baker partner names: hGH, hGHbp
#PubMed IDs: 9571026, 7504735
   {'PDBFileID' : '1A22',  'PartnerDefinitionID' : 1,
       'LName' : 'Human growth hormone ', 'LShortName' : 'hGH', 'LHTMLName' : 'hGH',
       'RName' : 'Human growth hormone binding protein', 'RShortName' : 'hGHbp', 'RHTMLName' : 'hGHbp',
       'IsDimerStructure' : None,
       'LChains': ['A'],
       'RChains': ['B'],
       },
         
#1A4Y
#Kortemme & Baker partner names: Angiogenin, Rnase Inh
   {'PDBFileID' : '1A4Y',  'PartnerDefinitionID' : 1,
       'LName' : 'RNase inhibitor', 'LShortName' : 'hRI', 'LHTMLName' : 'hRI',
       'RName' : 'Angiogenin', 'RShortName' : 'Ang', 'RHTMLName' : 'Ang',
       'IsDimerStructure' : True,
       'LChains': ['A'],
       'RChains': ['B'],
       },
        
#1AHW
#The PDB also has duplicate chains D, E, and F. Tanja used chains A (Light), B (Heavy), and C (Tissue Factor)
#Kortemme & Baker partner names: TF
   {'PDBFileID' : '1AHW',  'PartnerDefinitionID' : 1,
       'LName' : 'Tissue Factor', 'LShortName' : 'TF', 'LHTMLName' : 'TF',
       'RName' : 'Fab 5G9', 'RShortName' : '5G9', 'RHTMLName' : '5G9',
       'IsDimerStructure' : True,
       'LChains': ['C'],
       'RChains': ['A', 'B'],
       },
         
#1BRS
#The PDB also has duplicate chains B, E, C, and F
#Kortemme & Baker partner names: Barnase, Barstar
   {'PDBFileID' : '1BRS',  'PartnerDefinitionID' : 1,
       'LName' : 'Barnase', 'LShortName' : 'Barnase', 'LHTMLName' : 'Barnase',
       'RName' : 'Barstar', 'RShortName' : 'Barstar', 'RHTMLName' : 'Barstar',
       'IsDimerStructure' : True,
       'LChains': ['A'],
       'RChains': ['D'],
       },
         
#1BXI
#Kortemme & Baker partner names: Im9
#PubMed IDs: 9425068
   {'PDBFileID' : '1BXI',  'PartnerDefinitionID' : 1,
       'LName' : 'Colicin E9 DNase', 'LShortName' : 'E9 DNase', 'LHTMLName' : 'E9 DNase',
       'RName' : 'Colicin Immunity Protein Im9', 'RShortName' : 'Im9', 'RHTMLName' : 'Im9',
       'IsDimerStructure' : True,
       'LChains': ['B'],
       'RChains': ['A'],
       },
#1CBW
#Kortemme & Baker partner names: BPTI
   {'PDBFileID' : '1CBW',  'PartnerDefinitionID' : 1,
       'LName' : 'Bovine chymotrypsin', 'LShortName' : 'Chymotrypsin', 'LHTMLName' : 'Chymotrypsin',
       'RName' : 'Basic pancreatic trypsin inhibitor', 'RShortName' : 'BPTI', 'RHTMLName' : 'BPTI',
       'IsDimerStructure' : True,
       'LChains': ['A', 'B', 'C'],
       'RChains': ['D'],
       },
       
#1DAN
#PubMed IDs: 7654692
#            http://pubs.acs.org/doi/abs/10.1021/bi00033a009
#Kortemme & Baker partner names: TF
   {'PDBFileID' : '1DAN',  'PartnerDefinitionID' : 1,
       'LName' : 'Human tissue factor', 'LShortName' : 'TF', 'LHTMLName' : 'TF',
       'RName' : 'Human factor VIIa', 'RShortName' : 'FVIIa', 'RHTMLName' : 'FVIIa',
       'IsDimerStructure' : True,
       'LChains': ['T', 'U'],
       'RChains': ['L', 'H'],
       },
         
#1DFJ
#Kortemme & Baker partner names: Rnase Inh
   {'PDBFileID' : '1DFJ',  'PartnerDefinitionID' : 1,
       'LName' : 'Ribonuclease A', 'LShortName' : 'RNase A', 'LHTMLName' : 'RNase A',
       'RName' : 'Ribonuclease inhibitor', 'RShortName' : 'RI', 'RHTMLName' : 'RI',
       'IsDimerStructure' : True,
       'LChains': ['E'],
       'RChains': ['I'],
       },
         
#1DN2
#PubMed IDs: 10678837
#Kortemme & Baker partner names: IgG, Peptide
   {'PDBFileID' : '1DN2',  'PartnerDefinitionID' : 1,
       'LName' : 'Human immunoglobulin G', 'LShortName' : 'IgG', 'LHTMLName' : 'IgG',
       'RName' : 'Peptide', 'RShortName' : 'Peptide', 'RHTMLName' : 'Peptide',
       'IsDimerStructure' : None,
       'LChains': ['A'],
       'RChains': ['E'],
       },
         
#1F47
#PubMed IDs: 10880432
#Kortemme & Baker partner names: FTSZ fragm.
   {'PDBFileID' : '1F47',  'PartnerDefinitionID' : 1,
       'LName' : 'Z interacting protein A', 'LShortName' : 'ZipA', 'LHTMLName' : 'ZipA',
       'RName' : 'Filamenting temperature-sensitive mutant Z', 'RShortName' : 'FtsZ', 'RHTMLName' : 'FtsZ',
       'IsDimerStructure' : None,
       'LChains': ['B'],
       'RChains': ['A'],
       },
         
#1FC2
#Kortemme & Baker partner names: Protein A
# Use IgG Fc or IgG1 or hFc?
   {'PDBFileID' : '1FC2',  'PartnerDefinitionID' : 1,
       'LName' : 'Human Fc fragment of immunoglobulin G', 'LShortName' : ' hFc', 'LHTMLName' : 'hFc',
       'RName' : 'Fragment B of Protein A', 'RShortName' : 'Protein A', 'RHTMLName' : 'Protein A',
       'IsDimerStructure' : True,
       'LChains': ['D'],
       'RChains': ['C'],
       },
         
#1FCC
#PubMed IDs: 10452608
#Kortemme & Baker partner names: Protein G
   {'PDBFileID' : '1FCC',  'PartnerDefinitionID' : 1,
       'LName' : 'B1 domain (C2 fragment) of streptococcal protein G', 'LShortName' : 'B1', 'LHTMLName' : 'B1',
       'RName' : 'Human Fc fragment of immunoglobulin G', 'RShortName' : 'hFc', 'RHTMLName' : 'hFc',
       'IsDimerStructure' : None,
       'LChains': ['C'],
       'RChains': ['A'],
       },
         
#1GC1
#Kortemme & Baker partner names: CD4
   {'PDBFileID' : '1GC1',  'PartnerDefinitionID' : 1,
       'LName' : 'CD4 receptor', 'LShortName' : 'CD4', 'LHTMLName' : 'CD4',
       'RName' : 'HIV gp120 envelope glycoprotein', 'RShortName' : 'gp120', 'RHTMLName' : 'gp120',
       'IsDimerStructure' : True,
       'LChains': ['C'],
       'RChains': ['G'],
       },
        
#1JCK
#Kortemme & Baker partner names: SEC3
#ASEdb names: SEC3-1A4, TCR Vb
   {'PDBFileID' : '1JCK',  'PartnerDefinitionID' : 1,
       'LName' : 'T-cell receptor beta-chain', 'LShortName' : 'TCR', 'LHTMLName' : 'TCR',
       'RName' : 'Staphylococcus aureus enterotoxins C3 superantigen', 'RShortName' : 'SEC3', 'RHTMLName' : 'SEC3',
       'IsDimerStructure' : True,
       'LChains': ['A', 'C'],
       'RChains': ['B'],
       },
        
#1JRH
#PubMed IDs: 11123892, 9878445
#Kortemme & Baker partner names: A6, Interferon ?
   {'PDBFileID' : '1JRH',  'PartnerDefinitionID' : 1,
       'LName' : 'A6 fab-IFNgammaR1-108 complex', 'LShortName' : 'A6', 'LHTMLName' : 'A6',
      'RName' : 'Interferon gamma receptor (IFNgammaR) alpha-chain', 'RShortName' : 'hIFNgR', 'RHTMLName' : 'hIFNR',
       'IsDimerStructure' : None,
       'LChains': ['L', 'H'],
       'RChains': ['I'],
       },
         
#1NMB  FAB NC10 antibody HL, N9 influenza virus neuraminidase antigen N
#PubMed IDs: 9579662
#Kortemme & Baker partner names: NC10
   {'PDBFileID' : '1NMB',  'PartnerDefinitionID' : 1,
       'LName' : 'NC10 antibody', 'LShortName' : 'NC10 scFv', 'LHTMLName' : 'NC10 scFv',
       'RName' : 'Influenza virus glycoprotein neuraminidase', 'RShortName' : 'NA', 'RHTMLName' : 'NA',
       'IsDimerStructure' : None,
       'LChains': ['L', 'H'],
       'RChains': ['N'],
       },
         
#1VFB
#Kortemme & Baker partner names: D1.3, HEL
   {'PDBFileID' : '1VFB',  'PartnerDefinitionID' : 1,
       'LName' : 'Antibody D1.3', 'LShortName' : 'D1.3', 'LHTMLName' : 'D1.3',
       'RName' : 'Hen egg-white lysozyme', 'RShortName' : 'HEL', 'RHTMLName' : 'HEL',
       'IsDimerStructure' : True,
       'LChains': ['A', 'B'],
       'RChains': ['C'],
       },
         
#2PTC
#Kortemme & Baker partner names: BPTI
   {'PDBFileID' : '2PTC',  'PartnerDefinitionID' : 1,
       'LName' : 'Beta trypsin', 'LShortName' : 'Trypsin', 'LHTMLName' : 'Trypsin',
       'RName' : 'Bovine pancreatic trypsin inhibitor', 'RShortName' : 'BPTI', 'RHTMLName' : 'BPTI',
       'IsDimerStructure' : True,
       'LChains': ['E'],
       'RChains': ['I'],
       },
         
#3HFM
#PubMed IDs: 10338006
#Kortemme & Baker partner names: HEL, HYHEL-10
   {'PDBFileID' : '3HFM',  'PartnerDefinitionID' : 1,
       'RName' : 'Monoclonal antibody HyHEL-10', 'RShortName' : 'HyHEL-10', 'RHTMLName' : 'HyHEL-10',
       'LName' : 'Hen egg-white lysozyme', 'LShortName' : 'HEL', 'LHTMLName' : 'HEL',
       'IsDimerStructure' : True,
       'LChains': ['L', 'H'],
       'RChains': ['Y'],
       },
]
