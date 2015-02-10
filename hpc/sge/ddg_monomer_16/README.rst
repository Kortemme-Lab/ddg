==============================
Sun Grid Engine job submission
==============================

The protocols/ddg_monomer_16 script functions as both a submission script for a workstation and for an SGE cluster
engine. In the documentation below, we only describe the changes needed to run the script on the default dataset on an
SGE cluster; the user is referred to the documentation in protocols/ddg_monomer_16/README.rst for a more detailed
explanation of what the script does and how to change the input parameters.

In particular, the *Configuration* section of protocols/ddg_monomer_16/README.rst should be followed for the commands below
to work.

----------------------------
Paths and extensions
----------------------------

The command lines below use placeholders for paths and extensons. Please change these appropriately *e.g.*:

::

  WORKING_DIRECTORY=.
  BENCHMARK_PATH=<path/to/sequence-tolerance>
  OUTPUT_DIRECTORY=<directory created by the preminimization step>

The output directory will be named according the the current date and username *e.g.* 15-02-02-12-00_username_ddg_monomer_16.

-----------------------------
How to run the full benchmark
-----------------------------

These commands run the full benchmark for the different datasets.

**Kellogg dataset**

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python preminimization_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_ddg.py --force
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python ddg_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_analysis.py --force

**Guerois dataset**

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py -d ${BENCHMARK_PATH}/input/json/guerois.json
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python preminimization_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_ddg.py --force
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python ddg_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_analysis.py --force

**Potapov dataset**

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py -d ${BENCHMARK_PATH}/input/json/potapov.json
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python preminimization_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_ddg.py --force
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python ddg_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_analysis.py --force

**ProTherm\ \* dataset**

::

  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_preminimization.py -d ${BENCHMARK_PATH}/input/json/curatedprotherm.json
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python preminimization_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_ddg.py --force
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16/job_output/${OUTPUT_DIRECTORY}/
  python ddg_step.py

  # if execution is successful
  cd ${BENCHMARK_PATH}/protocols/ddg_monomer_16
  python run_analysis.py --force

