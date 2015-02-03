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

By default, the scripts are set up to run the benchmark using the Kellogg dataset.

----------------------
(Pre)minimization step
----------------------

The first step of the protocol is to run preminimization. The examples below do a test run which can be used to quickly
check that the installation was successful. To run a full benchmark, omit '--test' in the command lines below. To see
more options, run *python run_preminimization.py -h*.

::

  cd protocols/ddg_monomer_16
  python run_preminimization.py --test

This will create the job_output folder and a subfolder for the test run e.g. job_output/15-02-02-12-00_username_ddg_monomer_16.
The preminimization step is then run as follows:

::

  cd job_output/15-02-02-12-00_username_ddg_monomer_16/
  python cstmin_run.py

This creates preminimized structures used for the |DDG| step in the job_output/15-02-02-12-00_username_ddg_monomer_16/preminimization.

----------
|DDG| step
----------

The next step of the protocol is to run preminimization.

::

  cd protocols/ddg_monomer_16
  python run_preminimization.py --test

This will create the job_output folder and a subfolder for the test run e.g. job_output/15-02-02-12-00_username_ddg_monomer_16.
The preminimization step is then run as follows:

::

  cd job_output/15-02-02-12-00_username_ddg_monomer_16/
  python cstmin_run.py

This creates preminimized structures used for the |DDG| step in the job_output/15-02-02-12-00_username_ddg_monomer_16/preminimization.




.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G