====================================
Metrics
====================================

This benchmark uses three different metrics:

- Pearson's correlation coefficient;
- Mean absolute error (MAE); and
- Fraction correct.

The metrics are not mutually exclusive but each has a separate focus. Pearson's R indicates the level of correlation
between experimental and predicted values but ignores the average errors in cases. The MAE reports this error which
is important when using a DDG protocol for protein design. Finally, the fraction correct measures how likely we are to
correctly predict hotspots or neutral mutations.

For certain applications, the user may be more interested in one or two of the metrics above however the combination of
all three metrics reports how useful a DDG protocol is in general.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pearson's correlation coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pearson's correlation coefficient measures the linear correlation between a set of X values and their corresponding Y values.
In our analysis, we are measuring the correlation between experimentally determined |DDG| values for a set of mutations and
the computationally predicted values. The range of values is from -1.0 to 1.0. A good positive correlation is indicated by
a value between 0.7 and 1.0. A good negative correlation is indicated by a value between -0.7 and -1.0.  Of particular
note is that the Pearson correlation coefficient is invariant to the scale of the Y (predicted) values.

~~~~~~~~~~~~~~~~~~~~~~~~~
Mean absolute error (MAE)
~~~~~~~~~~~~~~~~~~~~~~~~~

The mean absolute error for this benchmark is defined as the mean of the absolute differences between experimental and predicted |DDG|
values. Note that the MAE can be large even if a correlation is high as its value is related to the scale of the predicted values.
Therefore, the MAE is an important metric if a |DDG| protocol is used in protein design as the individual predicted values
are of more importance in this context rather than the relative predictive accuracy (as measured by Pearson's coefficient).

If the scoring function is scaled such that one score unit is roughly equivalent to kcal/mol then a good value for MAE would
be a value within the range of experimental error.

~~~~~~~~~~~~~~~~
Fraction correct
~~~~~~~~~~~~~~~~

The fraction correct metric measures the fraction of dataset cases where a mutation was correctly classified as (de)stabilizing
or neutral. The metric is defined as the number of correctly classified mutations divided by the number of mutations in the
benchmark set so values range from 0.0 to 1.0.

The range of values for neutrality may vary; Kellogg et al. define neutrality as values between -1 kcal/mol
and 1 kcal/mol, exclusive, for experimental DDG values which we adopt here. They also define neutrality for predicted
values using their protocol 16, included here in protocols/ddg_monomer_16, as approximately between -3 and 1.15 Rosetta energy
units (see Supporting Information, Figure 4(d) [1]_). We define neutrality as between -1.0 and 1.0 energy units
in the included statistical analysis so the values of this metric may differ in this benchmark capture for the Rosetta
protocol.

Depending on the definition of neutrality, it is possible to get a relatively high value for this metric with a set of
random predicted values. Therefore this metric, while a useful metric for reporting whether a method can correctly classify
the stability of a mutant, should be considered alongside the correlation and MAE.

================
Running analysis
================

The analyze.py script in this directory performs basic analysis, reporting the values of the metrics above.

More detailed analysis can be performed by implementing Python classes which interact with the classes defined in the
reports directory. An example implementation is provided in protocols/ddg_monomer_16/run_analysis.py. These reporting
classes break the analysis down over different criteria and may be useful in development. They rely on a number of tools
including R_ and pandas_.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Paths and extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command lines below use placeholders for paths and extensons. Please change these appropriately *e.g.*:

::

  WORKING_DIRECTORY=.
  BENCHMARK_PATH=<path/to/sequence-tolerance>
  OUTPUT_DIRECTORY=<directory created by the preminimization step>

The output directory will be named according the the current date and username *e.g.* 15-02-02-12-00_username_ddg_monomer_16.

~~~~~~~~~~~~~~
Dependencies
~~~~~~~~~~~~~~

The analysis scripts require the `R software suite <http://www.r-project.org>`_. The scripts have been tested using R
versions 2.12.1 and 3.1.1. They also require the following Python libraries:

 - numpy
 - scipy

~~~~~~~~~~~~~
Main analysis
~~~~~~~~~~~~~

Analysis is performed by a Python script analysis/analyze.py which also uses R to create a scatterplot. The input to the script should be a JSON or a commas-/tabs-
separated file (CSV/TSV).

------------------
Input file formats
------------------

If the JSON format is used, the object should be a list of associative arrays/dicts each of which has both an Experimental and a
Predicted field with floating-point values *e.g.*:

::

  [{'Experimental': 0.8, 'ID': '71689', 'Predicted': 2.14764},
   {'Experimental': 2.6, 'ID': '71692', 'Predicted': 3.88848},
   ...
   {'Experimental': 1.4, 'ID': '76748', 'Predicted': 4.9911}]

The ID field above is unnecessary; any fields besides the required two fields will be ignored.

If the CSV/TSV format is used, all non-empty file lines not containing data should be preceded with a '#' character. The
first two columns will be used for the Experimental and Predicted values respectively *e.g.*:

::

  #Experimental,Predicted,ID,RecordID
  0.800000,2.147640,71689,332
  2.600000,3.888480,71692,333
  ...
  1.400000,4.991100,76748,944

----------------
Running analysis
----------------

::

  cd ${BENCHMARK_PATH}/analysis
  python analyze.py ${BENCHMARK_PATH}/output/sample/kellogg_r57471.txt

The analysis script prints out statistics on the input dataset, including the metrics above, and creates a scatterplot. By
default, this scatterplot will be named scatterplot.png in the current working directory but the --output argument can be used to specify the filename and file type (PNG or PDF format) *e.g.*

::

  cd ${BENCHMARK_PATH}/analysis
  python analyze.py ${BENCHMARK_PATH}/output/sample/kellogg_r57471.txt --output myplot.pdf  


-----------------
Rosetta protocols
-----------------


Note that it is not necessary to call the analyze.py explicitly for the included Rosetta protocols as their analysis scripts call the underlying statistical functions. See the relevant documentation in the protocols subdirectories for more details.

.. |Dgr|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |ring|  unicode:: U+002DA .. RING ABOVE
.. |DDGH2O| replace:: |Dgr|\ |Dgr|\ G H\ :sub:`2`\ O
.. |DDG| replace:: |Dgr|\ |Dgr|\ G


.. [1] Kellogg, EH, Leaver-Fay, A, Baker, D. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. 2011. Proteins. 79(3):830-8. `doi: 10.1002/prot.22921 <https://dx.doi.org/10.1002/prot.22921>`_.

.. _R: https://www.r-project.org

.. _pandas: http://pandas.pydata.org