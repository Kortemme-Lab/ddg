====================================
Input data
====================================

-------------
DDG benchmark
-------------

The input data for the DDG benchmark are DDG values determined from experimental assays and published in the literature.
The data are presented in the form of curated datasets published in the literature for the purposes of benchmarking DDG protocols.
Most of the data in the datasets appears to be derived from the ProTherm database [1]_, [2]_ which contains a large amount
of thermodynamic parameters from published experimental assays. ProTherm also stores the references associated with
these parameters which greatly enriches the data.

There are four curated datasets in this benchmark capture, compiled by separate teams. The datasets overlap but
differ in their choice of DDG values, methods of denaturation (*chemical, thermal, *etc.*), and structural (PDB) files.
Some of these choices depend on when the dataset was published, particularly: i) the choice of PDB file as higher-resolution
structures have since been deposited in the PDB; and ii) the set of published DDG values. We refer the
reader to the original publications for the full details of the selection criteria.

Note: The three previously published datasets are not presented in their original forms here. In particular, we have made
changes to the list of mutations, the PDB ID, or the DDG values in cases where we suspected that the published values may
be incorrect or in cases where the PDB IDs have since been deprecated. Our changes are hopefully but not necessarily correct
and we welcome any corrections from the community.

[discuss mapping of datasets to ProTherm records]

--------------------------
Alanine scanning benchmark
--------------------------

- AlaScan-GPK
- SKEMPI

-----------------
Benchmark formats
-----------------

The dataset files are provided in `JavaScript Object Notation <http://www.json.org/>`_ (JSON) and/or comma-separated values
(CSV) format. When both formats are available, the JSON files contain more detailed information whereas we have tried to
keep the CSV files as terse as possible.

We chose to use both formats for flexibility. JSON files can be easily incorporated into applications written in, but not
limited to, C, C++, C#, Java, JavaScript, Perl, and Python. They are human-readable but verbose. CSV files may be more
immediately read in a text editor or spreadsheet application but may require more work to incorporate into an application.

==============
PDB structures
==============


- input/pdbs
- input/hydrogen_pdbs
- pdbs.json
- pdbs.csv


==========
References
==========

- references.json
- references.csv

========
Datasets
========

The following datasets are presented herein:


---------------------------
Guerois et al. [3]
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/guerois.json, input/csv/guerois.csv

---------------------------
Potapov et al. [4]
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/potapov.json, input/csv/potapov.csv

---------------------------
Kellogg et al. [4]
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/kellogg.json, input/csv/kellogg.csv

---------------------------
ProTherm*
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/curatedprotherm.json, input/csv/curatedprotherm.csv

---------------------------
AlaScan-GPK
---------------------------

Records: ...
Unique PDB IDs: ...
Files: input/json/alascan-gpk.json, input/csv/alascan-gpk.csv

---------------------------
References
---------------------------

.. [1] Gromiha, MM, An, J, Kono, H, Oobatake, M, Uedaira, H, Sarai, A. ProTherm: Thermodynamic Database for Proteins and Mutants. 1999. Nucl. Acids Res. 27(1):286-288. `doi: 10.1093/nar/27.1.286 <http://dx.doi.org/10.1093/nar/27.1.286>`_.

.. [2] Kumar, SK, Bava, KA, Gromiha, MM, Prabakaran, P, Kitajima, K, Uedaira, H, Sarai, A. ProTherm and ProNIT: thermodynamic databases for proteins and protein–nucleic acid interactions. 2006. Nucleic Acids Res. 34(Database issue):D204-6. `doi: 10.1093/nar/gkj103 <http://dx.doi.org/10.1093/nar/gkj103>`_.
