================
Alanine scanning
================

-------------------
General Information
-------------------

Created by: Tanja Kortemme, David Baker [1]_

Software suite: Rosetta

Protocol directory: protocols/alanine-scanning

-----------
Description
-----------

The previously published alanine scanning protocol (see `References`_ below) has been recreated in modern Rosetta using Rosetta scripts.
Unlike the generalized |DDG| protocol, which minimizes the entire protein structure, the alanine scanning protocol avoids any perturbation of the backbone or side chains, other than the residue being mutated, which is placed into a low-energy rotamer using the Rosetta “packer”.
This minimal perturbation relies on the fact that the overall protein structure is unlikely to change much after a single point mutation to alanine, making the input crystal structure a good approximation for the mutant structure.

As the alanine scanning protocol does not perturb the protein backbone or side chains (other than the mutant residue), this protocol is not suitable for use on mutations outside of the interface.
A mutation outside of the interface will result in a negligible change in total score without the use of a more intensive sampling protocol.

As in |DDG|, the metrics used to measure success in this benchmark are: i) the linear correlation (Pearson coefficient) between experimental and predicted values; ii) the mean absolute error (MAE) of same; and iii) the stability classification accuracy, which measures whether a mutation was correctly predicted to be stabilizing, destabilizing, or neutral.

-------------
Configuration
-------------

The file protocols/rosetta/settings.json should be created and configured according to the user's system. An
example configuration file in provided in protocols/rosetta/settings.json.example. The *local_rosetta_installation_path*
setting should be set to the path of a Rosetta installation on a workstation. The *cluster_rosetta_installation_path* path
is used when submitting a job on a computational cluster (scripts are included for running the protocol on an SGE cluster - see
hpc/sge/alanine-scanning). *rosetta_binary_type* should be set to the suffix used in the Rosetta application binaries [2]_.

The two locations of the Rosetta installations should be the directories which contains the checkout of the Rosetta repository
*i.e.* the directories containing:

- database
- README.md
- source
- tests

The file uses the
`JSON <http://www.json.org/>`_ format.

---------------------------------
Included benchmark helper scripts
---------------------------------

There is a helper script included in this folder which can be used to run the benchmark with minimal setup.

**setup_alanine_scanning.py**

When run with no arguments, this script sets up a benchmarking run on the full alanine scanning set with the standard protocol.

The helper setup script generates a job output directory containing the run script *alascan_run.py*.
This script can be directly run (using all available local processors) or submitted to a SGE cluster via the qsub command.
As the protocol runs very quickly, the entire benchmark can be run on an average workstation on the scale of minutes.

Expected execution time per alanine point mutation: 5-30 seconds

setup_alanine_scanning.py needs to be configured to run correctly on a user's system. This is explained in `Configuration`_.

----------------------------
Paths and extensions
----------------------------

The command lines below use placeholders for paths and extensons. Please change these appropriately *e.g.*:

::

  WORKING_DIRECTORY=.
  BENCHMARK_PATH=<path/to/protocols/alanine-scanning>
  OUTPUT_DIRECTORY=<directory created by the preminimization step>

The output directory will be named according the the current date and username *e.g.* *15-02-02-12-00_username_alanine_scanning*.

------------
Dependencies
------------

The benchmark and analysis scripts use `Python <https://www.python.org/>`_ and the `R software suite <http://www.r-project.org>`_ respectively. These
scripts have been tested with:

- Python 2.6.6
- Python 2.7.8 and R 3.1.1

The analysis scripts also require the following Python libraries:

- numpy
- scipy

And the folllowing R library:

- ggplot2

An installation of Rosetta is required for this method. Rosetta can be downloaded `here <https://www.rosettacommons.org/>`__
and is freely available for academic use. Details of how to install Rosetta can be found in the `User Guide <https://www.rosettacommons.org/docs/latest/>`__.


=====================
Running the benchmark
=====================

We first run a setup script creates the input directories, the output directories, and an execution/run script.
This execution script is then called to run that step.

The examples below run the benchmark on the entire benchmark dataset, with default options.

----------
Setup step
----------

The first step of the protocol generates a self-contained output directory with all files needed to run the benchmark.

::

  cd ${BENCHMARK_PATH}/protocols/alanine-scanning
  python setup_alanine_scanning.py

This will create the default folder, *job_output*, and a subfolder for the test run *e.g.* job_output/${OUTPUT_DIRECTORY}.

--------------------
Run (execution) step
--------------------

Alanine scanning is then run as follows:

::

  cd ${BENCHMARK_PATH}/protocols/alanine-scanning/job_output/${OUTPUT_DIRECTORY}/
  python alascan_run.py

This step completes the protocol and outputs each subtask within the bencmark into a new subdirectory with the following structure:

::

   pdb_id/score_function_name/pdb_id_resnum

The |DDG| values of mutation to alanine are saved in multiple places in each subdirectory:

- At the bottom of the output PDB file
- As part of the output in *rosetta.out.gz* (grep for lines containing protocols.features.DdGFeatures)
- Within the Sqlite3 *rosetta_output.db3* files (in the ddG table)

We can now run the analysis script to complete the benchmark run.

--------
Analysis
--------

The analyze.py script is used to parse the outputted |DDG| values from the outputted Sqlite3 databases and produce statistics and plots.
The script takes one argument: the output directory contained the data to be processed:
::

  cd ${BENCHMARK_PATH}/protocols/alanine-scanning
  python analyze.py ${BENCHMARK_PATH}/protocols/alanine-scanning/job_output/${OUTPUT_DIRECTORY}

This script creates three kinds of files in the output directory:

- PDF files plotting the experimental alanine scanning |DDG| values vs. predicted values.
  Each PDF file is named according to which score function was used to generate the predictions
  (the interface weights, score12 score function, and talaris score function are used by default).
- run_name-stats.txt contains the benchmark metrics that are also printed to standard output (example below)
- A csv file containing the experimental data and all predictions, for later analysis

Example benchmark statistics output (saved in run_name-stats.txt and printed to console):

::

   interface_pts - ddg_obs vs tanja_ddg_calc
   Fraction correct                : 0.742
   Fraction correct (fuzzy)        : 0.747
   Gamma correlation coef.         : 0.388
   Kolmogorov-Smirnov test (XY)    : 0.112 (2-tailed p-value=0.102418547943)
   MAE                             : 1.054
   Pearson's R                     : 0.516 (2-tailed p-value=2.94319824392e-17)
   Spearman's R                    : 0.547 (2-tailed p-value=1.33174136685e-19)
   X-axis Kolmogorov-Smirnov test  : 0.401 (p-value=0.0)
   X-axis normality test           : 58.487 (2-sided chi^2 p-value=1.99371969776e-13)
   Y-axis Kolmogorov-Smirnov test  : 0.470 (p-value=0.0)
   Y-axis normality test           : 126.341 (2-sided chi^2 p-value=3.67548452526e-28)

   interface_pts - ddg_obs vs interface
   Fraction correct                : 0.687
   Fraction correct (fuzzy)        : 0.698
   Gamma correlation coef.         : 0.359
   Kolmogorov-Smirnov test (XY)    : 0.159 (2-tailed p-value=0.00486842162526)
   MAE                             : 1.059
   Pearson's R                     : 0.488 (2-tailed p-value=2.59582916372e-15)
   Spearman's R                    : 0.513 (2-tailed p-value=5.18644180117e-17)
   X-axis Kolmogorov-Smirnov test  : 0.401 (p-value=0.0)
   X-axis normality test           : 58.487 (2-sided chi^2 p-value=1.99371969776e-13)
   Y-axis Kolmogorov-Smirnov test  : 0.418 (p-value=0.0)
   Y-axis normality test           : 53.026 (2-sided chi^2 p-value=3.0591076075e-12)

The first section of statistics shows the results for the original Rosetta implementation of alanine scanning for comparison.
The second (and each following) show the results for each score function used; results for interface weighting is shown in this example.

================================
Appendix A: Command line options
================================

The setup script above has additional command-line options for specifying the location of the output files, among other things.
The help files for these options can be viewed by using the -h or --help flag e.g.

::

  cd ${BENCHMARK_PATH}/protocols/alanine-scanning
  python setup_alanine_scanning.py --help
  python setup_alanine_scanning.py -h

The help text is generated automatically from the Python scripts using the `docopt <https://github.com/docopt>`_ module.

==========
References
==========

Kortemme, T, Baker, D. A simple physical model for binding energy hot spots in protein–protein complexes.
Proc Natl Acad Sci U S A. 2002 Oct 29;99(22):14116-21. Epub 2002 Oct 15.
`doi: 10.1073/pnas.202485799 <https://dx.doi.org/10.1073/pnas.202485799>`_.

Kortemme T, Kim DE, Baker D. Computational alanine scanning of protein-protein interfaces.
Sci STKE. 2004 Feb 3;2004(219):pl2.
`doi: 10.1126/stke.2192004pl2 <https://dx.doi.org/10.1126/stke.2192004pl2>`_.

.. [1] The Rosetta application was written by the authors above. This protocol capture was compiled by Kyle Barlow. Any errors in the protocol capture are likely to be the fault of the compiler rather than that of the original authors. Please contact support@kortemmelab.ucsf.edu with any issues which may arise.
.. [2] By default, a Linux release build of Rosetta built with GCC will append the suffix '.linuxgccrelease' to binaries *e.g.* ddg_monomer.linuxgccrelease is the binary for the backrub application.

.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G
