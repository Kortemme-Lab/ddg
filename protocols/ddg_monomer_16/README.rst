=================================================
Protein stability protocol 1: ddg_monomer, row 16
=================================================

-------------------
General Information
-------------------

Created by: Elizabeth Kellogg, Andrew Leaver-Fay, David Baker [1]_

Software suite: Rosetta

Protocol directory: protocols/ddg_monomer_16

-----------
Description
-----------

The ddg_monomer Rosetta application was created and benchmarked using twenty different sampling strategies ((see `References`_
below). The protocol listed as number 16 in that paper performed well in general and is used as the basis for this protocol
capture, hence *ddg_monomer_16*. Protocol 16 combines a soft-repulsive potential for conformational sampling of side-chains with a standard
hard-repulsive potential for minimization to achieve higher prediction accuracy, following the observation that predictive
methods are more accurate when the resolution of the force field is matched to the granularity of the sampling method.
This approach helps to prevent sampling from being driven into local minima in the energy landscape [ref Dantas, Kellogg].

There are two main computational steps in the protocol - a preminimization step and a |DDG| step. The steps use the Rosetta
minimize_with_cst and ddg_monomer applications respectively, both of which are described `here <https://www.rosettacommons.org/docs/latest/ddg-monomer.html>`__

In preminimization, the input structure is stripped to a single monomer chain and minimized. This creates two types of output -
the minimized protein structure and a set of constraints used to prevent the backbone from moving too much in the next step.

In the DDG step, the minimized structures and constraints are fed into the ddg_monomer application which creates fifty
pairs of wildtype and mutant structures using the sampling strategy described above, calculates their score using the
Rosetta full-atom force field, and reports the difference (predicted |DDG|) as the difference between the best-scoring
wild-type structure and the best-scoring mutant structure.

The output from the DDG step is then processed by the analysis scripts which compare the predicted |DDG| values against
published |DDG| values from the literature, generating: i) a correlation coefficient; ii) a mean absolute error; and iii) the fraction
of cases ([0.0, 1.0]) that were predicted correctly *i.e.* whether a mutation was correctly predicted to be (de)stabilizing
or neutral.

----------------
Included scripts
----------------

There are two scripts included in this folder which can be used to help run the benchmark.

**run_preminimization.py**

Expected execution time: 1-2 minutes per structure

This script runs the preminimization step of the protocol. The inputs to the protocol are the dataset of experimental DDG
values and an output directory, which default to the Kellogg dataset and the protocols/ddg_monomer_16/job_output directory
respectively.

This script determines the set of PDB chains used in the dataset and then creates a preminimizated, monomeric structure
for each of these chains. The output of the preminimization execution is stored in rosetta.out.gz for use in the run_ddg.py
script.

run_preminimization.py needs to be configured to run correctly on a user's system. This is explained in `Configuration`_.

**run_ddg.py**

Expected execution time: <30 hours per dataset record if the length of the protein is <400 residues. [2]_

This script runs the main computational step of the protocol. First, a set of constraints is generated from the output
created by the preminimization step. Next, fifty pairs of wild-type and mutant structures are generated and the |DDG| value
is determined as the difference between the best-scoring (using the Rosetta scoring function) wild-type structure and the
best-scoring mutant structure.


-------------
Configuration
-------------

The file protocols/rosetta/settings.json should be created and configured according to the user's system. An
example configuration file in provided in protocols/rosetta/settings.json.example. The *local_rosetta_installation_path*
setting should be set to the path of a Rosetta installation on a workstation. The *cluster_rosetta_installation_path* path
is used when submitting a job on a computational cluster (scripts are included for running the protocol on an SGE cluster - see
hpc/sge/ddg_monomer_16). *rosetta_binary_type* should be set to the suffix used in the Rosetta application binaries [3]_.

The two locations of the Rosetta installations should be the directories which contains the checkout of the Rosetta repository
*i.e.* the directories containing:

- database
- README.md
- source
- tests

The file uses the
`JSON <http://www.json.org/>`_ format.


----------------------------
Paths and extensions
----------------------------

The command lines below use placeholders for paths and extensons. Please change these appropriately *e.g.*:

::

  WORKING_DIRECTORY=.
  BENCHMARK_PATH=<path/to/sequence-tolerance>
  OUTPUT_DIRECTORY=<directory created by the preminimization step>

The output directory will be named according the the current date and username *e.g.* 15-02-02-12-00_username_ddg_monomer_16.

------------
Dependencies
------------

The benchmark and analysis scripts use `Python <https://www.python.org/>`_ and the `R software suite <http://www.r-project.org>`_ respectively. These
scripts have been tested with:
- Python 2.6.6
- Python 2.7.8 and R 3.1.1

An installation of Rosetta is required for this method. Rosetta can be downloaded `here <https://www.rosettacommons.org/>`__
and is freely available for academic use. Details of how to install Rosetta can be found in the `User Guide <https://www.rosettacommons.org/docs/latest/>`__.


=====================
Running the benchmark
=====================

In both steps, we first run a setup script creates the input directories, the output directories, and an execution script.
This execution script is then called to run that step.

The examples below run the benchmark using the default dataset, the Kellogg set. This can be overridden by using the
-d flag in the *preminimization* step to provide another suitable JSON input file *e.g.* '-d ${BENCHMARK_PATH}/input/json/potapov.json'.


----------------------
(Pre)minimization step
----------------------

The first step of the protocol generates preminimized monomeric structures and sets of constraints for all of the protein
chains in the dataset.

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py

This will create the default folder, *job_output*, and a subfolder for the test run *e.g.* job_output/${OUTPUT_DIRECTORY}.
The preminimization step is then run as follows:

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python preminimization_step.py

This creates preminimized structures used for the |DDG| step in the job_output/${OUTPUT_DIRECTORY}/preminimization. A
copy of the dataset JSON file is stored in job_output/${OUTPUT_DIRECTORY}/ for use in the following
steps.

As mentioned above, the benchmarking dataset is chosen at this stage of execution and defaults to the Kellogg dataset. The
various |DDG| datasets can be set up to run as follows:

::

  # Benchmark the Kellogg dataset
  python run_preminimization.py
  # Benchmark the Guerois dataset
  python run_preminimization.py -d ${BENCHMARK_PATH}/input/json/guerois.json
  # Benchmark the Potapov dataset
  python run_preminimization.py -d ${BENCHMARK_PATH}/input/json/potapov.json
  # Benchmark the ProTherm* dataset
  python run_preminimization.py -d ${BENCHMARK_PATH}/input/json/curatedprotherm.json

----------
|DDG| step
----------

The next step of the protocol is to run ddg_monomer. If preminimization was run in the default output folder (job_output) then
the run_ddg.py script prompts the user to ask whether the most recent subfolder should be used. This prompt can be skipped
by using the --force argument (as used below). If preminimization was run in a different folder, this should be supplied to the
script via the -o option.

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_ddg.py --force

This sets up the input files for the run in the same directory as used in the preminimization step. The |DDG| step is then run as follows:

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python ddg_step.py

This step completes the protocol and outputs pairs (50 pairs by default) of wildtype and mutant structures and |DDG| scores for
each record in the input dataset. Each record has a RecordID field in the dataset JSON file. The output for the dataset record
with RecordID n is stored in the directory ddg/n.

We can now run the analysis script to complete the benchmark run.

--------
Analysis
--------

Before we can run the analysis script analsis/analyze.py, we need to compile the results of the benchmark run. This is
done with the run_analysis.py script, which also invokes analyze.py for convenience. As in the last step, if the default
output folder (job_output) was used for the first two steps then the run_analysis.py script prompts the user to ask
whether the most recent subfolder should be used. Again, this prompt can be skipped by using the --force argument. If
the benchmark was run in a different folder, this should be supplied to the script via the -o option.

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_analysis.py --force

This script creates four files in the output directory:

- analysis_input.json, a JSON file which contains experimental and predicted |DDG| values and dataset record IDs (to help identify outliers). This is then passed to analysis/analyze.py;
- analysis_input.csv, a CSV version of analysis_input.json;
- benchmark_data.json, a JSON file containing all of the Rosetta score components for the wildtype and mutant structures generated by the |DDG| step of the protocol. This is provided for convenience in case users wish to perform their own analysis;
- scatterplot.png [4]_, a scatterplot image plotting the experimental and predicted |DDG| values.

The analysis script also prints out the benchmark metrics to the terminal as well as a number of other metrics which may
also be of interest e.g.

::

  ********** Statistics **********
  Fraction correct                : 0.400
  Fraction correct (fuzzy)        : 0.405
  Gamma correlation coef.         : 0.326
  Kolmogorov-Smirnov test (XY)    : 0.150 (2-tailed p-value=0.965484740899)
  MAE                             : 1.703
  Pearson's R                     : 0.426 (2-tailed p-value=0.0613856027581)
  Spearman's R                    : 0.499 (2-tailed p-value=0.025021864609)
  X-axis Kolmogorov-Smirnov test  : 0.433 (p-value=0.000627758702818)
  X-axis normality test           : 0.003 (2-sided chi^2 p-value=0.998403992043)
  Y-axis Kolmogorov-Smirnov test  : 0.364 (p-value=0.00699618500741)
  Y-axis normality test           : 0.246 (2-sided chi^2 p-value=0.884305721297)

(Note that these statistics were generated from a test run - see below).

---------
Test mode
---------


Before starting a full benchmark run, it is advisable to make sure that everything works by running a test version of the
benchmark. In the test benchmark (for the Kellogg dataset), three preminimized structures are created corresponding to 20
records in the dataset. For each record, only 2 pairs of wildtype and mutant structures are generated. As such, results
produced in test mode should be ignored.

Test mode is enabled by passing the --test flag to both the run_preminimization.py and the run_ddg.py scripts. For example,
the command lines for an entire test run are as follows:

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py --test
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python preminimization_step.py

  [if execution is successful]
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_ddg.py --force --test
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python ddg_step.py

  [if execution is successful]
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_analysis.py --force


(To run the full benchmark, omit --test in the commands above)

================================
Appendix A: Command line options
================================

The scripts above have additional command-line options for specifying the location of the output files. The help files for
these options can be viewed by using the -h or --help flag e.g.

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py --help
  python run_ddg.py -h
  python run_analysis.py --help

For convenience, the options are printed below however we suggest that the --help flag is used in case this documentation
is not updated with changes to the code.

The help text is generated automatically from the Python scripts using the `docopt <https://github.com/docopt>`_ module.

----------------------
(Pre)minimization step
----------------------

Usage:
    run_preminimization.py [options]...

Options:

    -d --dataset DATASET
        A filepath to the input dataset in JSON format. [default: ../../input/json/kellogg.json]

    -o --output_directory OUTPUT_DIR
        The path where output data will be created. Output will be created inside a time-stamped subfolder of this directory. [default: ./job_output]

    --run_identifier RUN_ID
        A suffix used to name the output directory.

    --test
        When this option is set, a shorter version of the benchmark will run with fewer input structures, less fewer DDG experiments, and fewer generated structures. This should be used to test the scripts but not for analysis.

----------
|DDG| step
----------

Usage:
    run_ddg.py [options]...

Options:

    -o --output_directory OUTPUT_DIR
        The path to a directory previously generated from the run_preminimization script. This defaults to the most recent directory in job_output, if this exists.

    -n --num_struct NUM_STRUCT
        This specifies the number of wildtype/mutant structures generated. If this is used with --test then the --test value for this option takes priority. [default: 50]

    --force
        When this option is set, the most recent directory in job_output, if it exists, will be used without prompting the user.

    --test
        When this option is set, a shorter version of the benchmark will run with fewer input structures, less fewer DDG experiments, and fewer generated structures. This should be used to test the scripts but not for analysis.

----------------------
Analysis
----------------------

Usage:
    run_analysis.py [options]...

Options:

    -o --output_directory OUTPUT_DIR
        The path to a directory previously generated from the run_preminimization script. This defaults to the most recent directory in job_output, if this exists.

    -p --scatterplot_filename SCATTERPLOT_FILE
        The filename of the scatterplot to be generated in the output directory (unless --skip_analysis is set). [default: scatterplot.png]

    --force
        When this option is set, the most recent directory in job_output, if it exists, will be used without prompting the user.

    --skip_analysis
        When this option is set, the analysis script is not invoked once the analysis files are generated.


==========
References
==========

Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. 2011. Proteins. 79(3):830-8. `doi: 10.1002/prot.22921 <https://dx.doi.org/10.1002/prot.22921>`_.


.. [1] The Rosetta application was written by the authors above. This protocol capture was compiled by Shane O'Connor. Any errors in the protocol capture are likely to be the fault of the compiler rather than that of the original authors. Please contact support@kortemmelab.ucsf.edu with any issues which may arise.
.. [2] Given the amount of computational resources needed for benchmarking using this protocol, we recommend that this step is performed using cluster, grid, or cloud computing. The execution time is proportional to the number of wildtype/mutant pairs generated which is 50 by default. This number can be reduced but we would recommend using at least the default value.
.. [3] By default, a Linux release build of Rosetta built with GCC will append the suffix '.linuxgccrelease' to binaries *e.g.* ddg_monomer.linuxgccrelease is the binary for the backrub application.
.. [4] This default filename can be overridden using the --scatterplot_filename option of the run_analysis.py script.




.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G