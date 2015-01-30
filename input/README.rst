====================================
Input data
====================================

---------------
|DDG| benchmark
---------------

The input data for the |DDG| benchmark are |DDG| values determined from experimental assays and published in the literature.
The data are presented in the form of curated datasets published for the purposes of benchmarking |DDG| protocols.

The datasets are divided into two classes:
 - monomeric datasets; each data point specifies a single protein chain, a set of mutations on that chain, and an experimental |DDG| value measuring protein stability;
 - protein-protein interface datasets; each data point specifies the two partners of the interface (where each partner may consist of multiple protein chains), a set of mutations on the chains, and an experimental |DDG| value measuring binding affinity.

Each |DDG| value in the datasets is presented with references to published material describing the conditions in which the
measurement was taken. We provide a bibliography containing details of the references.

------------
File formats
------------

The dataset files, reference file, and PDB description files are provided in `JavaScript Object Notation <http://www.json.org/>`_ (JSON)
and/or comma-separated values (CSV) format. When both formats are available, the JSON files contain more detailed information
whereas we have tried to keep the CSV files as terse as possible.

We chose to use both formats for flexibility. JSON files can be easily incorporated into applications written in, but not
limited to, C, C++, C#, Java, JavaScript, Perl, and Python. They are human-readable and well-suited for web applications
but are verbose. CSV files may be more immediately read in a text editor or spreadsheet application but may require additional
programming to incorporate into an application.

====================================
Monomeric datasets (thermostability)
====================================

There are five curated monomeric datasets in this benchmark capture, compiled by separate teams. The datasets overlap but
differ in their choice of |DDG| values, methods of denaturation (chemical, thermal, *etc.*), and structural (PDB) files.
Some of these choices depend on when the dataset was published, particularly: i) the choice of PDB file as higher-resolution
structures have since been deposited in the PDB; and ii) the set of published |DDG| values. We refer the
reader to the original publications for the full details of the selection criteria.

The three previous published datasets use data compiled in the ProTherm database [1]_, [2]_ which contains a large amount of
thermodynamic parameters from published experimental assays. ProTherm also stores the references associated with
these parameters which greatly enriches the data. The datasets were originally published using |DDG| values with no explicit
connection to the original sources in the literature but we have made efforts to add this data back into the datasets by
mapping the dataset records between themselves and the ProTherm database.

ProTherm distinguishes between types of |DDG| measurements - in general |DDG| and |DDGH2O| values respectively denote thermal
and denaturant denaturation.

Note: The three previously published datasets are not presented in their original forms here. In particular, we have made
changes to the list of mutations, the PDB ID, or the |DDG| values in cases where we suspected that the published values may
be incorrect or in cases where the PDB IDs have since been deprecated. Our changes are hopefully - but not necessarily -
correct. We welcome any corrections from the community.

---------------------------------
Guerois et al. [3]_ dataset, 2002
---------------------------------

Records: 1005

Unique PDB IDs: 81

Files: input/json/guerois.json, input/csv/guerois.csv

This dataset was compiled for benchmarking the FOLDEF application. ProTherm was a source for many of the datapoints. We use
the combined dataset (training dataset and blind test dataset) of monomeric single mutants used for testing protein stability
prediction. The dataset strikes a balance in terms of the secondary structure type of the wildtype residue,
hydrophobic-hydrophobic mutations versus changes in polarity, and buried mutations versus exposed mutations.

Each record appears to correspond with one recorded |DDG| value in ProTherm.

---------------------------
Potapov et al. [4]_ dataset
---------------------------

Records: 2154

Unique PDB IDs: 84

Files: input/json/potapov.json, input/csv/potapov.csv

This dataset was compiled for benchmarking different methods for measuring protein stability. The data originates from
single mutations in both the Guerois dataset, which contains many ProTherm records, and the ProTherm database directly.

Most records appears to correspond with one recorded |DDG| value in ProTherm but approximately 19% of the records appear
to use a mean value from multiple recorded |DDG| values.

As with all of the previous published datasets, we have modified this dataset. In particular, 10 of the records in the
dataset presented herein have duplicate records with the same PDB ID and mutation. These are due to our replacement of
1BKS with 1WQ5, 2LZMA with 1L63, and 1HGU with 3HHR.

---------------------------
Kellogg et al. [5]_ dataset
---------------------------

Records: 1210

Unique PDB IDs: 75

Files: input/json/kellogg.json, input/csv/kellogg.csv

This dataset was compiled for benchmarking 20 configurations of the Rosetta ddg_monomer application. The data is taken
from single mutations in the ProTherm database although some records seem to be taken from the Guerois dataset or else added independently.
Higher-resolution PDB structures were used when a choice existed and structures with more than 350 residues were omitted
to reduce the computational cost of the benchmark. The lowest-resolution PDB structure in the dataset is 2.8A.

Most records appears to correspond with one recorded |DDG| value in ProTherm but approximately 6% of the records appear
to use a mean value from multiple recorded |DDG| values.

-----------------
ProTherm* dataset
-----------------

Records: 2971

Unique PDB IDs: 119

Files: input/json/curatedprotherm.json, input/csv/curatedprotherm.csv

This dataset was compiled for this benchmark capture. It is a curated subset of single mutations in the ProTherm dataset
with the intention of using as many ProTherm records as possible. This dataset is therefore biased according to the results
published in the literature. The selection criteria are:

- any mutations in transmembrane protein chains are omitted;
- |DDGH2O| values are favored over |DDG| values when available;
- all PDB structures are determined or partially determined by X-ray defraction with a resolution of at least 2.5A;
- records where any two |DDG| values vary by more than 2.5 kcal/mol are omitted.

Approximately 28% of the records use a mean value from multiple recorded |DDG| values.

-----------
AlaScan-GPK
-----------

Records: 768

Unique PDB IDs: 56

Files: input/json/alascan-gpk.json, input/csv/alascan-gpk.csv

This dataset consists of all of the point mutations to alanine that are present in the Guerois, Potapov, and Kellogg datasets.
The method of construction was as follows:

- consider all point mutations to alanine in the union of the datasets;
- for each mutation, take:

 - the set of |DDGH2O| values from ProTherm if available, otherwise the set of |DDG| values. If this set contained |DDG| values used in the datasets then we took the mean value of the intersection otherwise we took the mean value of the entire set;
 - PDB structures determined by X-ray defraction over those determined by NMR, if available;
 - the highest resolution PDB structure used in the datasets.

Thus, the |DDG| values and PDB IDs may differ from the original datasets for some records. Approximately 12% of the records use a mean value from multiple recorded |DDG| values.


=====================================================
Protein-protein interface datasets (binding affinity)
=====================================================

todo: describe SKEMPI

-----------
SKEMPI [6]_
-----------

Records: todo:

Unique PDB IDs: todo:

Files: todo:, todo:


==============
PDB structures
==============

The PDB structures used for the benchmark are provided in input/pdbs.

todo: describe input/hydrogen_pdbs

The resolutions and methods of determination (*e.g.* X-ray defraction, NMR, *etc.*) for the PDB structures are listed
in input/json/pdbs.json and input/csv/pdbs.csv.

=======================
Dataset reference files
=======================

The datasets tie each experimental |DDG| value to a reference using an identifying string (typically a PubMed ID). The details
of these references - authors, title, publication, publication date - are provided in input/json/references.json and
input/csv/references.csv. The references for the novel datasets presented herein - the ProTherm* and AlaScan-GPK
datasets - should be accurate. The references for the Guerois, Potapov, and Kellogg datasets have been deduced so they may
not be entirely accurate.

==========
References
==========

.. [1] Gromiha, MM, An, J, Kono, H, Oobatake, M, Uedaira, H, Sarai, A. ProTherm: Thermodynamic Database for Proteins and Mutants. 1999. Nucl. Acids Res. 27(1):286-288. `doi: 10.1093/nar/27.1.286 <https://dx.doi.org/10.1093/nar/27.1.286>`_.

.. [2] Kumar, SK, Bava, KA, Gromiha, MM, Prabakaran, P, Kitajima, K, Uedaira, H, Sarai, A. ProTherm and ProNIT: thermodynamic databases for proteins and protein–nucleic acid interactions. 2006. Nucleic Acids Res. 34(Database issue):D204-6. `doi: 10.1093/nar/gkj103 <https://dx.doi.org/10.1093/nar/gkj103>`_.

.. [3] Guerois, R, Nielsen, JE, Serrano, L. Predicting changes in the stability of proteins and protein complexes: a study of more than 1000 mutations. 2002. J Mol Biol. 320(2):369-87. `doi: 10.1016/S0022-2836(02)00442-4 <https://dx.doi.org/10.1016/S0022-2836(02)00442-4>`_.

.. [4] Potapov, V, Cohen, M, Schreiber, G. Assessing computational methods for predicting protein stability upon mutation: good on average but not in the details. 2009. Protein Eng Des Sel. 22(9):553-60. `doi: 10.1093/protein/gzp030 <https://dx.doi.org/10.1093/protein/gzp030>`_.

.. [5] Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. 2011. Proteins. 79(3):830-8. `doi: 10.1002/prot.22921 <https://dx.doi.org/10.1002/prot.22921>`_.

.. [6] Moal, IH, Fernández-Recio, J. SKEMPI: a Structural Kinetic and Energetic database of Mutant Protein Interactions and its use in empirical models. 2012. Bioinformatics. 28(20):2600-7. `doi: 10.1093/bioinformatics/bts489 <https://dx.doi.org/10.1093/bioinformatics/bts489>`_.

.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G
