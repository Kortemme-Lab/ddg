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

---------------------------
Guerois et al. [3]_ dataset
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/guerois.json, input/csv/guerois.csv

---------------------------
Potapov et al. [4]_ dataset
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/potapov.json, input/csv/potapov.csv

---------------------------
Kellogg et al. [5]_ dataset
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/kellogg.json, input/csv/kellogg.csv

-----------------
ProTherm* dataset
-----------------

Records: ...
Unique PDB IDs: ...
Files: input/json/curatedprotherm.json, input/csv/curatedprotherm.csv

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
 - the highest resolution PDB structure used in the datasets.

Thus, the |DDG| values and PDB IDs may differ from the original datasets for some records.

=====================================================
Protein-protein interface datasets (binding affinity)
=====================================================

todo: describe SKEMPI

-----------
SKEMPI [6]_
-----------

Records: ...
Unique PDB IDs: ...
Files: todo:, todo:



==============
PDB structures
==============

The PDB structures used for the benchmark are provided in input/pdbs.
todo: describe input/hydrogen_pdbs

The resolutions and methods of determination (*e.g.* X-ray crystallograhy, NMR, *etc.*) for the PDB structures are listed
in the input/json/pdbs.json and input/csv/pdbs.csv.

=======================
Dataset reference files
=======================

The datasets tie each experimental |DDG| value to a reference using an identifying string (typically a PubMed ID). The details
of these references - authors, title, publication, publication date - are provided in input/json/references.json and
input/csv/references.csv


==========
References
==========

.. [1] Gromiha, MM, An, J, Kono, H, Oobatake, M, Uedaira, H, Sarai, A. ProTherm: Thermodynamic Database for Proteins and Mutants. 1999. Nucl. Acids Res. 27(1):286-288. `doi: 10.1093/nar/27.1.286 <http://dx.doi.org/10.1093/nar/27.1.286>`_.

.. [2] Kumar, SK, Bava, KA, Gromiha, MM, Prabakaran, P, Kitajima, K, Uedaira, H, Sarai, A. ProTherm and ProNIT: thermodynamic databases for proteins and protein–nucleic acid interactions. 2006. Nucleic Acids Res. 34(Database issue):D204-6. `doi: 10.1093/nar/gkj103 <http://dx.doi.org/10.1093/nar/gkj103>`_.

.. [3] Guerois

.. [4] Potapov

.. [5] Kellogg

.. [6] SKEMPI

.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G
