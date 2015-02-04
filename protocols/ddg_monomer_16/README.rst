=====================
Running the benchmark
=====================

The ddg_monomer_16 protocol is run in two steps. In the first step, the input structures are pruned to the protein chain
of interest and are then minimized (using the minimize_with_cst Rosetta application). This creates two types of output -
the minimized protein structures and a set of constraints. These outputs are used in the second step, the |DDG| computation.
In the second step, the minimized structures and constraints are fed into the ddg_monomer application which creates a set
pairs of wildtype and mutant structures, calculates their score using the Rosetta full-atom force field, and reports the
difference (predicted |DDG|) between these pair of scores. The final |DDG| score is given as the difference in the lowest
energy pair of structures.

The examples below do a test run which can be used to quickly check that the installation was successful. To run a full
benchmark, omit '--test' in the command lines below. To see the options for the scripts, run:

::

  python run_preminimization.py -h
  python run_ddg.py -h

In both steps, we first run a script to set up the directories and create an execution script. We then run that execution
script.

Note: The scripts are set up to run the benchmark using the Kellogg dataset by default. This can be overridden by using the
-d flag in the preminimization step to provide another suitable JSON input file *e.g.* '-d ../../input/json/potapov.json'.

----------------------
(Pre)minimization step
----------------------

The first step of the protocol is to generate preminimized monomeric structures and sets of constraints.

::

  cd protocols/ddg_monomer_16
  python run_preminimization.py --test

This will create the default folder, *job_output*, and a subfolder for the test run *e.g.* job_output/15-02-02-12-00_username_ddg_monomer_16.
The preminimization step is then run as follows:

::

  cd job_output/15-02-02-12-00_username_ddg_monomer_16/
  python preminimization_step.py

This creates preminimized structures used for the |DDG| step in the job_output/15-02-02-12-00_username_ddg_monomer_16/preminimization. A
copy of the dataset JSON file is stored in job_output/15-02-02-12-00_username_ddg_monomer_16/ for use in the following
steps.

----------
|DDG| step
----------

The next step of the protocol is to run ddg_monomer. If preminimization was run in the default output folder (job_output) then
the run_ddg.py script prompts the user to ask whether the most recent subfolder should be used. This prompt can be skipped
by using the --force argument (as used below). If preminimization was run in a different folder, this should be supplied to the
script via the -o option.

::

  cd protocols/ddg_monomer_16
  python run_ddg.py --force --test

This sets up the input files for the run.

::

  cd job_output/15-02-02-12-00_username_ddg_monomer_16/
  python ddg_step.py

This step completes the protocol and outputs pairs (50 by default) of wildtype and mutant structures and |DDG| scores for
each record in the input dataset. The output for the record with id record_id is stored in the directory ddg/record_id.
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

  cd protocols/ddg_monomer_16
  python run_analysis.dy --force

This script creates benchmark_output.csv and benchmark_output.json in the output directory.

The next step of the protocol is to run ddg_monomer. If preminimization was run in the default folder (job_output) then
the run_ddg.py script prompts the user to ask whether the most recent subfolder should be used. This prompt can be skipped
by using the -f argument (as used below). If preminimization was run in a different folder, this should be supplied to the
script via the -o option.

::

  cd protocols/ddg_monomer_16
  python run_ddg.py -f --test


.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G