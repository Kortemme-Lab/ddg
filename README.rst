========================================
Energetic effects of mutation benchmarks
========================================

This benchmark capture contains two related benchmarks, both of which measure the accuracy with which a protocol predicts
the energetic effects of mutation.

-----------------------------------
|DDG| (protein stability) benchmark
-----------------------------------

Monomeric |DDG| protocols predict the change in protein stability which occurs as a result of mutagenesis. This benchmark
includes four curated datasets of experimentally measured values which can be used to evaluate the accuracy of a protocol.

This benchmark includes:

- four curated datasets, three of which were previously published (see input/README.rst). Most of the data was taken from the `ProTherm database <http://www.abren.net/protherm>`_;
- scripts to run a Rosetta protocol described in Kellogg et al. (2011);
- an analysis script that output the metrics used for analysis. The script also outputs a scatterplot plotting experimental |DDG| values in kcal/mol against predicted values (in whichever scoring unit is used by the protocol);
- sample output data which can be used to test the analysis scripts.

Benchmarking a complete dataset is quite computationally intensive so we recommend that a benchmark run be performed on a cluster, grid, or cloud computing resource. We have provided scripts to run this benchmark on a Sun Grid Engine cluster (see hpc/sge/ddg_monomer_16/README.rst).

--------------------------
Alanine scanning benchmark
--------------------------

This benchmark tests the ability of computational alanine scanning protocols to recapitulate the results of experimental
alanine scanning, wherein each residue in a range of interest is individually mutated to the “neutral” residue alanine,
disrupting any hydrogen bonding or shape-specific interactions made by the wild-type side chain without major disruptions
to secondary structure. This technique allows a direct before and after interrogation of the energetics of folding, binding,
or other functional effects of a single residue position through comparison of the properties of each alanine mutant to those
of the wild-type protein *i.e.* alanine scanning can be used to ascertain the contribution of each individual residue to
these interaction energies.

This benchmark includes:

- a previously published set of 233 mutations (to alanine) in 19 different protein-protein interfaces with known crystal structures (see `Kortemme & Baker, 2002 <References>`_);
- scripts to run a new RosettaScripts protocol which has been designed to emulate the protocol described in Kortemme & Baker (2002);
- an analysis script that output the metrics used for analysis. The script also outputs a scatterplot plotting experimental |DDG| values against predicted values (in whichever scoring unit is used by the protocol);

The RosettaScripts alanine scanning protocol is not computationally intensive so this benchmark can be performed on a typical laptop or workstation.

---------
Licensing
---------

The contents of the repository *where possible* are licensed under the MIT License. The license only applies to files which
either: i) include the license statement; or ii) which are explicitly listed in some file in the repository as being covered
by the license. All other files may be covered under a separate license. The LICENSE file in the root of this repository
is present only for the convenience of the user to indicate the license which covers any novel content presented herein.

-------------------------
Downloading the benchmark
-------------------------

The benchmark is hosted on GitHub. The most recent version can be checked out using the `git <http://git-scm.com/>`_ command-line tool:

::

  git clone https://github.com/Kortemme-Lab/ddg.git

---------------------------
Directories in this archive
---------------------------

This archive contains the following directories:

- *input* : contains the input files for the benchmarks. These take the form of PDB files and datasets of experimental |DDG| values. The input files are described in more detail in input/README.rst;
- *output/sample* : contains sample output data that can be used to test the stand-alone analysis script;
- *analysis* : contains the stand-alone analysis script (analyze.py) and the analysis functions (stats.py). All protocols are expected to produce output that will work with both scripts. These scripts do not need to be called directly - each benchmark has a separate analysis step which performs analysis;
- *protocols* : contains the scripts needed to run the benchmarks. The scripts for each protocol are provided in a specific subdirectory;
- *protocols/alanine-scanning* : contains the scripts needed to run the alanine scanning benchmark using a Rosetta protocol;
- *protocols/ddg_monomer_16* : contains the scripts needed to run the protein stability benchmark using the Rosetta ddg_monomer protocol;
- *hpc* : contains scripts that can be used to run the entire benchmark using specific cluster architectures. For practical reasons, a limited number of cluster systems are supported. Please feel free to provide scripts which run the benchmark for your particular cluster system.
- *hpc/sge/ddg_monomer_16* : contains scripts that can be used to run the the protein stability benchmark using ddg_monomer on a Sun Grid Engine cluster.


=========
Protocols
=========

This repository contains a protocol which can be used to run the |DDG| benchmark and another which can be used to run the
alanine scanning benchmark. We welcome the inclusion of more protocols. Please contact support@kortemmelab.ucsf if you wish
to contribute towards the repository.

Each protocol is accompanied by specific documentation in its protocol directory.

-------------------------------------------------
Protein stability protocol 1: ddg_monomer, row 16
-------------------------------------------------

Created by: Elizabeth Kellogg, Andrew Leaver-Fay, David Baker [1]_

Software suite: Rosetta

Protocol directory: protocols/ddg_monomer_16

----------------------------------------------------
Alanine scanning protocol 1: RosettaScripts protocol
----------------------------------------------------

Created by: Kyle Barlow

Software suite: Rosetta

Protocol directory: protocols/alanine-scanning


========
Analysis
========

The same set of analysis scripts is used by all protocols. Conceptually, the analysis scripts should be a black box that
is separated from the output of each protocol by an interface. The expected input format is described in analysis/README.rst.

The analysis scripts generates three metrics which can be used to evaluate the results of the |DDG| and alanine scanning
simulations and also produces a scatterplot of the experimental and predicted values. The benchmark analysis is described
in more detail in analysis/README.rst.


==========
References
==========

-----------------------------------
|DDG| (protein stability) benchmark
-----------------------------------

Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. 2011.
Proteins. 79(3):830-8. `doi: 10.1002/prot.22921 <https://dx.doi.org/10.1002/prot.22921>`_.

--------------------------
Alanine scanning benchmark
--------------------------

Kortemme, T, Baker, D. A simple physical model for binding energy hot spots in protein–protein complexes.
Proc Natl Acad Sci U S A. 2002 Oct 29;99(22):14116-21. Epub 2002 Oct 15.
`doi: 10.1073/pnas.202485799 <https://dx.doi.org/10.1073/pnas.202485799>`_.

Kortemme T, Kim DE, Baker D. Computational alanine scanning of protein-protein interfaces.
Sci STKE. 2004 Feb 3;2004(219):pl2.
`doi: 10.1126/stke.2192004pl2 <https://dx.doi.org/10.1126/stke.2192004pl2>`_.


=====
Notes
=====

.. [1] The Rosetta application was written by the authors above. This protocol capture was compiled by Shane O'Connor. Any errors in the protocol capture are likely to be the fault of the compiler rather than that of the original authors. Please contact support@kortemmelab.ucsf.edu with any issues which may arise.


.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G


